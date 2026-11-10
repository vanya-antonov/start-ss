#!/usr/bin/env python3
"""Retrieve human orthologs in the nine Ensembl 116 mammalian GTF species."""

import argparse
import csv
import json
import logging
import os
import tempfile
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data" / "mammals"
DEFAULT_INPUT = DATA_DIR / "01_Homo_sapiens.canonical_transcripts.tsv"
DEFAULT_OUTPUT = DATA_DIR / "03_human_orthologs.tsv"
CACHE_DIR = DATA_DIR / ".cache" / "03_orthologs"
ENSEMBL_REST = "https://rest.ensembl.org"
HUMAN_SPECIES = "homo_sapiens"
REQUEST_INTERVAL = 0.1  # 10 requests/second, below the Ensembl REST limit.
MAX_ATTEMPTS = 5
TIMEOUT_SECONDS = 60
LOG = logging.getLogger("identify_orthologs")

# These are the nine nonhuman species with GTFs in the documented Ensembl 116 set.
TARGET_SPECIES = {
    "mouse": "mus_musculus",
    "cow": "bos_taurus",
    "sheep": "ovis_aries",
    "panda": "ailuropoda_melanoleuca",
    "elephant": "loxodonta_africana",
    "dog": "canis_lupus_familiaris",
    "horse": "equus_caballus",
    "rabbit": "oryctolagus_cuniculus",
    "marmoset": "callithrix_jacchus",
}
ORTHOLOGY_TYPES = {
    "ortholog_one2one": 0,
    "ortholog_one2many": 1,
    "ortholog_many2many": 2,
}
HUMAN_COLUMNS = ("gene_id", "gene_name")


class EnsemblRequestError(RuntimeError):
    """Raised when an Ensembl REST request fails after retries."""


def stable_id(identifier):
    """Remove an optional Ensembl version suffix for REST stable-ID queries."""
    return identifier.strip().rsplit(".", maxsplit=1)[0]


def read_human_genes(input_path):
    """Read unique human gene IDs and names from the canonical-transcript TSV."""
    if not input_path.is_file():
        raise FileNotFoundError(f"Human transcript table not found: {input_path}")

    genes = {}
    with input_path.open(encoding="utf-8-sig", newline="") as input_file:
        reader = csv.DictReader(input_file, delimiter="\t")
        missing = set(HUMAN_COLUMNS) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                "Input TSV is missing required columns: " + ", ".join(sorted(missing))
            )
        for line_number, row in enumerate(reader, start=2):
            gene_id = (row.get("gene_id") or "").strip()
            if not gene_id:
                raise ValueError(f"Missing gene_id at input line {line_number}")
            query_id = stable_id(gene_id)
            gene_name = (row.get("gene_name") or "").strip()
            genes.setdefault(query_id, {"gene_id": query_id, "gene_name": gene_name})
            if not genes[query_id]["gene_name"] and gene_name:
                genes[query_id]["gene_name"] = gene_name
    if not genes:
        raise ValueError("Input TSV contains no human genes")
    return genes


def read_cache(cache_path):
    """Load a cached JSON response, returning None when it is not available."""
    try:
        with cache_path.open(encoding="utf-8") as cache_file:
            return json.load(cache_file)
    except FileNotFoundError:
        return None
    except (OSError, json.JSONDecodeError) as error:
        LOG.warning("Ignoring unreadable cache %s: %s", cache_path, error)
        return None


def write_json_atomically(path, content):
    """Write cached JSON through a temporary file to allow safe restarts."""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as output_file:
            json.dump(content, output_file, separators=(",", ":"))
        Path(temporary_name).replace(path)
    except Exception:
        Path(temporary_name).unlink(missing_ok=True)
        raise


class EnsemblClient:
    """Rate-limited Ensembl REST client with retry and response caching."""

    def __init__(self, cache_dir):
        self.cache_dir = cache_dir
        self.next_request_at = 0.0

    def get_json(self, endpoint, cache_name):
        cache_path = self.cache_dir / f"{cache_name}.json"
        cached = read_cache(cache_path)
        if cached is not None:
            return cached

        url = f"{ENSEMBL_REST}{endpoint}"
        for attempt in range(1, MAX_ATTEMPTS + 1):
            delay = self.next_request_at - time.monotonic()
            if delay > 0:
                time.sleep(delay)
            self.next_request_at = time.monotonic() + REQUEST_INTERVAL
            request = Request(
                url,
                headers={
                    "Accept": "application/json",
                    "User-Agent": "start-ss mammalian orthology analysis",
                },
            )
            retry_after = None
            try:
                with urlopen(request, timeout=TIMEOUT_SECONDS) as response:
                    data = json.loads(response.read().decode("utf-8"))
                write_json_atomically(cache_path, data)
                return data
            except HTTPError as error:
                retry_after = error.headers.get("Retry-After")
                retryable = error.code == 429 or 500 <= error.code <= 599
                message = f"HTTP {error.code}: {error.reason}"
            except (URLError, TimeoutError, json.JSONDecodeError) as error:
                retryable = True
                message = str(error)

            if not retryable or attempt == MAX_ATTEMPTS:
                raise EnsemblRequestError(
                    f"Ensembl request failed after {attempt} attempt(s): {url}: {message}"
                )
            try:
                wait_seconds = max(0.0, float(retry_after)) if retry_after else 0.0
            except ValueError:
                wait_seconds = 0.0
            if wait_seconds <= 0:
                wait_seconds = min(2 ** (attempt - 1), 30)
            LOG.warning(
                "Request attempt %d/%d failed (%s); retrying in %.1f seconds",
                attempt,
                MAX_ATTEMPTS,
                message,
                wait_seconds,
            )
            time.sleep(wait_seconds)
        raise AssertionError("Unreachable retry loop")

    def homologies_for_gene(self, gene_id):
        """Return cached or fetched Ensembl Compara orthology data for a gene."""
        endpoint = f"/homology/id/{HUMAN_SPECIES}/{gene_id}?{urlencode({'type': 'orthologues', 'sequence': 'none'})}"
        cache_name = f"homology_{gene_id}"
        return self.get_json(endpoint, cache_name)

    def gene_name(self, gene_id):
        """Look up a target gene's display name, if Ensembl provides one."""
        endpoint = f"/lookup/id/{gene_id}?{urlencode({'expand': 0})}"
        result = self.get_json(endpoint, f"lookup_{gene_id}")
        return (result.get("display_name") or "").strip()


def target_species_key(homology):
    """Return the short project prefix for a homology target species."""
    species = (homology.get("target") or {}).get("species", "").lower()
    for prefix, ensembl_species in TARGET_SPECIES.items():
        if species == ensembl_species:
            return prefix
    return None


def optional_number(*values):
    """Format the first available numeric coverage/identity value."""
    for value in values:
        if value is not None and value != "":
            return str(value)
    return ""


def homology_fields(homology):
    """Read Compara fields while tolerating fields absent from REST responses."""
    query = homology.get("query") or {}
    target = homology.get("target") or {}
    return {
        "gene_id": target.get("id", ""),
        "gene_name": target.get("display_name") or target.get("name") or "",
        "orthology_type": homology.get("type", ""),
        "perc_identity": optional_number(target.get("perc_id"), homology.get("perc_id")),
        "query_coverage": optional_number(
            query.get("perc_cov"), query.get("qcov"), homology.get("query_coverage")
        ),
        "target_coverage": optional_number(
            target.get("perc_cov"), target.get("tcov"), homology.get("target_coverage")
        ),
    }


def collect_orthologs(client, genes):
    """Fetch orthologies and target names, preserving every target gene."""
    collected = {}
    total = len(genes)
    for gene_number, (human_id, human) in enumerate(genes.items(), start=1):
        response = client.homologies_for_gene(human_id)
        data = response.get("data") or []
        by_species = {species: [] for species in TARGET_SPECIES}
        homologies = data[0].get("homologies", []) if data else []
        for homology in homologies:
            species = target_species_key(homology)
            if species is None:
                continue
            fields = homology_fields(homology)
            if not fields["gene_id"]:
                LOG.warning("Homology for %s in %s had no target ID", human_id, species)
                continue
            fields["gene_id"] = stable_id(fields["gene_id"])
            if not fields["gene_name"]:
                try:
                    fields["gene_name"] = client.gene_name(fields["gene_id"])
                except EnsemblRequestError as error:
                    LOG.warning("Could not retrieve name for %s: %s", fields["gene_id"], error)
            by_species[species].append(fields)

        for orthologs in by_species.values():
            orthologs.sort(
                key=lambda item: (
                    ORTHOLOGY_TYPES.get(item["orthology_type"], 99),
                    item["gene_id"],
                )
            )
        collected[human_id] = {"human": human, "orthologs": by_species}
        if gene_number % 100 == 0 or gene_number == total:
            LOG.info("Processed %d/%d human genes", gene_number, total)
    return collected


def output_columns():
    """Build the required output schema in stable species order."""
    columns = ["human_gene_name", "human_gene_id"]
    for species in TARGET_SPECIES:
        columns.extend(
            [
                f"{species}_gene_id",
                f"{species}_gene_name",
                f"{species}_orthology_type",
                f"{species}_perc_identity",
                f"{species}_query_coverage",
                f"{species}_target_coverage",
            ]
        )
    return columns


def output_rows(collected):
    """Build one output row per unique human stable gene ID."""
    rows = []
    for human_id in sorted(collected):
        entry = collected[human_id]
        row = {
            "human_gene_name": entry["human"]["gene_name"],
            "human_gene_id": human_id,
        }
        for species, orthologs in entry["orthologs"].items():
            for field in (
                "gene_id",
                "gene_name",
                "orthology_type",
                "perc_identity",
                "query_coverage",
                "target_coverage",
            ):
                row[f"{species}_{field}"] = ";".join(
                    ortholog[field] for ortholog in orthologs
                )
        rows.append(row)
    return rows


def write_output(output_path, rows):
    """Write the ortholog table atomically as a tab-separated file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output_path.name}.", dir=output_path.parent
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as output_file:
            writer = csv.DictWriter(
                output_file,
                fieldnames=output_columns(),
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)
        Path(temporary_name).replace(output_path)
    except Exception:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def summarize(collected):
    """Report per-species counts by orthology type and no detected ortholog."""
    total = len(collected)
    LOG.info("Total unique human genes processed: %d", total)
    for species in TARGET_SPECIES:
        counts = {orthology_type: 0 for orthology_type in ORTHOLOGY_TYPES}
        missing = 0
        for entry in collected.values():
            types = {item["orthology_type"] for item in entry["orthologs"][species]}
            if not types:
                missing += 1
            for orthology_type in counts:
                if orthology_type in types:
                    counts[orthology_type] += 1
        one_to_one = counts["ortholog_one2one"]
        percentage = (100 * one_to_one / total) if total else 0
        LOG.info(
            "%s: one-to-one %d (%.2f%%); one-to-many %d; many-to-many %d; "
            "without detected ortholog %d",
            species,
            one_to_one,
            percentage,
            counts["ortholog_one2many"],
            counts["ortholog_many2many"],
            missing,
        )


def main():
    parser = argparse.ArgumentParser(
        description="Identify human orthologous genes in nine mammalian species."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Human canonical transcript TSV")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Output ortholog TSV")
    parser.add_argument(
        "--cache-dir", type=Path, default=CACHE_DIR, help="Directory for restartable API response cache"
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    try:
        genes = read_human_genes(args.input)
        LOG.info("Read %d unique human genes from %s", len(genes), args.input)
        client = EnsemblClient(args.cache_dir)
        collected = collect_orthologs(client, genes)
        rows = output_rows(collected)
        write_output(args.output, rows)
    except (EnsemblRequestError, FileNotFoundError, OSError, ValueError) as error:
        parser.error(str(error))

    summarize(collected)
    LOG.info("Wrote %d rows to %s", len(rows), args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
