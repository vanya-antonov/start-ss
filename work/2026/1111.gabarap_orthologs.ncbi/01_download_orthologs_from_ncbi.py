#!/usr/bin/env python3
"""Download an NCBI Datasets ortholog set for a human gene symbol."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import sys
import time
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen


EXPERIMENT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = EXPERIMENT_DIR.parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "mammals"))
from constants import SPECIES_DICT  # noqa: E402


API_ROOT = "https://api.ncbi.nlm.nih.gov/datasets/v2"
HUMAN_TAXON_ID = SPECIES_DICT["human"]["ncbi_taxon_id"]
REQUEST_INTERVAL_SECONDS = 0.21  # 5 requests/s, NCBI's default limit.
MAX_ATTEMPTS = 5
TIMEOUT_SECONDS = 60
LOG = logging.getLogger("ncbi_gabarap_orthologs")
OUTPUT_COLUMNS = [
    "human_gene_symbol",
    "human_gene_id",
    "species",
    "taxon_id",
    "ortholog_gene_id",
    "ortholog_gene_symbol",
    "ortholog_type",
    "refseq_accession",
]


class NCBIRequestError(RuntimeError):
    """Raised when an NCBI Datasets request fails after retries."""


class NCBIClient:
    """Small rate-limited HTTP client that caches raw JSON responses."""

    def __init__(self, data_dir: Path) -> None:
        self.data_dir = data_dir
        self.next_request_at = 0.0

    def get_json(self, endpoint: str, cache_name: str) -> dict[str, Any]:
        """Get JSON from Datasets, retaining the exact decoded response."""
        path = self.data_dir / f"{cache_name}.json"
        if path.exists():
            with path.open(encoding="utf-8") as raw_file:
                return json.load(raw_file)

        url = f"{API_ROOT}{endpoint}"
        for attempt in range(1, MAX_ATTEMPTS + 1):
            delay = self.next_request_at - time.monotonic()
            if delay > 0:
                time.sleep(delay)
            self.next_request_at = time.monotonic() + REQUEST_INTERVAL_SECONDS
            request = Request(
                url,
                headers={
                    "Accept": "application/json",
                    "User-Agent": "start-ss NCBI orthology research workflow",
                },
            )
            try:
                with urlopen(request, timeout=TIMEOUT_SECONDS) as response:
                    payload = json.loads(response.read().decode("utf-8"))
                self.data_dir.mkdir(parents=True, exist_ok=True)
                with path.open("w", encoding="utf-8") as raw_file:
                    json.dump(payload, raw_file, indent=2, ensure_ascii=False)
                    raw_file.write("\n")
                return payload
            except HTTPError as error:
                retry_after = error.headers.get("Retry-After")
                retryable = error.code == 429 or 500 <= error.code < 600
                message = f"HTTP {error.code}: {error.reason}"
            except (URLError, TimeoutError, json.JSONDecodeError) as error:
                retry_after = None
                retryable = True
                message = str(error)
            if not retryable or attempt == MAX_ATTEMPTS:
                raise NCBIRequestError(f"NCBI request failed: {url}: {message}")
            try:
                wait_seconds = float(retry_after) if retry_after else 0.0
            except ValueError:
                wait_seconds = 0.0
            if wait_seconds <= 0:
                wait_seconds = min(2 ** (attempt - 1), 30)
            LOG.warning("Attempt %d failed (%s); retrying in %.1f s", attempt, message, wait_seconds)
            time.sleep(wait_seconds)
        raise AssertionError("Unreachable retry loop")


def configured_targets() -> list[tuple[str, str, int]]:
    """Return configured nonhuman species and taxonomy IDs in project order."""
    targets = []
    ordered_species = sorted(
        SPECIES_DICT.items(), key=lambda item: item[1]["sort_order"]
    )
    for key, details in ordered_species:
        if key == "human":
            continue
        taxon_id = details.get("ncbi_taxon_id")
        if not taxon_id:
            raise ValueError(f"No NCBI taxonomy ID configured for species key {key!r}")
        targets.append((key, details["species"], taxon_id))
    return targets


def reports_from_response(payload: dict[str, Any]) -> list[dict[str, Any]]:
    """Extract gene reports from documented Datasets response wrappers."""
    if "reports" in payload:
        reports = payload["reports"]
    elif "gene_reports" in payload:
        reports = payload["gene_reports"]
    else:
        raise ValueError("NCBI response did not contain a recognized gene reports array")
    if not isinstance(reports, list):
        raise ValueError("NCBI gene reports field was not a list")
    return [item.get("gene", item) for item in reports if isinstance(item, dict)]


def gene_fields(report: dict[str, Any]) -> dict[str, str]:
    """Normalize fields used by NCBI Gene reports."""
    gene_id = report.get("gene_id") or report.get("geneid") or ""
    symbol = report.get("symbol") or report.get("gene_symbol") or ""
    tax = report.get("tax_id") or report.get("taxid") or ""
    return {
        "gene_id": str(gene_id),
        "symbol": str(symbol),
        "tax_id": str(tax),
        "refseq_accession": refseq_accession(report),
        "ortholog_type": str(report.get("ortholog_type") or report.get("type") or ""),
    }


def refseq_accession(report: dict[str, Any]) -> str:
    """Choose an available RefSeq accession from common gene-report fields."""
    for key in ("refseq_accession", "refseq", "rna_accession"):
        value = report.get(key)
        if isinstance(value, str):
            return value
        if isinstance(value, list):
            for item in value:
                if isinstance(item, str) and item.startswith(("NM_", "NR_", "XM_", "XR_")):
                    return item
                if isinstance(item, dict):
                    accession = item.get("accession_version") or item.get("accession")
                    if accession:
                        return str(accession)
    products = report.get("transcripts") or report.get("gene_products") or []
    for product in products:
        if isinstance(product, dict):
            accession = product.get("accession_version") or product.get("accession")
            if accession and str(accession).startswith(("NM_", "NR_", "XM_", "XR_")):
                return str(accession)
    return ""


def resolve_human_gene(client: NCBIClient, symbol: str) -> dict[str, str]:
    """Resolve a symbol against human Gene and require one exact symbol match."""
    encoded = quote(symbol, safe="")
    endpoint = f"/gene/symbol/{encoded}/taxon/{HUMAN_TAXON_ID}/dataset_report"
    payload = client.get_json(endpoint, f"human_symbol_{safe_filename(symbol)}")
    candidates = [gene_fields(row) for row in reports_from_response(payload)]
    exact = [gene for gene in candidates if gene["symbol"].casefold() == symbol.casefold()]
    if len(exact) != 1:
        if not candidates:
            raise ValueError(f"No human NCBI Gene record found for symbol {symbol!r}")
        if not exact:
            raise ValueError(f"No exact human gene-symbol match for {symbol!r}; matches: {candidates}")
        raise ValueError(f"Ambiguous human gene symbol {symbol!r}: {exact}")
    if not exact[0]["gene_id"]:
        raise ValueError(f"Resolved record for {symbol!r} had no NCBI Gene ID")
    return exact[0]


def collect_ortholog_rows(
    client: NCBIClient, human_gene: dict[str, str], symbol: str
) -> list[dict[str, str]]:
    """Retrieve orthologs and emit at least one row per configured species."""
    endpoint = f"/gene/id/{quote(human_gene['gene_id'], safe='')}/orthologs"
    payload = client.get_json(endpoint, f"orthologs_{human_gene['gene_id']}")
    reports = reports_from_response(payload)
    targets = {str(taxon_id): (label, taxon_id) for _, label, taxon_id in configured_targets()}
    rows_by_taxon: dict[str, list[dict[str, str]]] = {taxon: [] for taxon in targets}
    for report in reports:
        fields = gene_fields(report)
        if fields["tax_id"] not in targets:
            continue
        species, taxon_id = targets[fields["tax_id"]]
        if not fields["gene_id"]:
            LOG.warning("Skipping ortholog report without GeneID for taxon %s", taxon_id)
            continue
        rows_by_taxon[fields["tax_id"]].append({
            "human_gene_symbol": symbol,
            "human_gene_id": human_gene["gene_id"],
            "species": species,
            "taxon_id": str(taxon_id),
            "ortholog_gene_id": fields["gene_id"] or "-",
            "ortholog_gene_symbol": fields["symbol"] or "-",
            "ortholog_type": fields["ortholog_type"] or "-",
            "refseq_accession": fields["refseq_accession"] or "-",
        })

    rows = []
    for taxon, (species, taxon_id) in targets.items():
        species_rows = rows_by_taxon[taxon]
        if species_rows:
            rows.extend(species_rows)
        else:
            rows.append({
                "human_gene_symbol": symbol,
                "human_gene_id": human_gene["gene_id"],
                "species": species,
                "taxon_id": str(taxon_id),
                "ortholog_gene_id": "-",
                "ortholog_gene_symbol": "-",
                "ortholog_type": "-",
                "refseq_accession": "-",
            })
    return rows


def safe_filename(symbol: str) -> str:
    """Make a symbol safe for use in an output filename."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", symbol.strip()).strip("._-") or "gene"


def write_table(path: Path, rows: list[dict[str, str]]) -> None:
    """Write ortholog rows to a tab-delimited file, including header if empty."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    """Run symbol resolution, ortholog retrieval, and table creation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gene_symbol", help="Human gene symbol, for example GABARAP")
    parser.add_argument(
        "output_file",
        nargs="?",
        help=(
            "optional output TSV path; relative paths are resolved from the "
            "current working directory"
        ),
    )
    args = parser.parse_args()
    symbol = args.gene_symbol.strip()
    if not symbol:
        parser.error("gene symbol cannot be empty")

    data_dir = EXPERIMENT_DIR / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    client = NCBIClient(data_dir)
    try:
        human_gene = resolve_human_gene(client, symbol)
        rows = collect_ortholog_rows(client, human_gene, symbol)
        if args.output_file:
            output_path = Path(args.output_file).expanduser()
            if not output_path.is_absolute():
                output_path = Path.cwd() / output_path
        else:
            results_dir = EXPERIMENT_DIR / "results"
            output_path = results_dir / f"01_{safe_filename(symbol)}_ncbi_orthologs.tsv"
        write_table(output_path, rows)
    except (NCBIRequestError, ValueError, OSError, json.JSONDecodeError) as error:
        LOG.error("Analysis stopped: %s", error)
        return 1

    print(f"Human {symbol} GeneID: {human_gene['gene_id']}")
    found_orthologs = sum(
        1 for row in rows
        if row["ortholog_gene_id"] != "-"
    )
    configured_species = len(configured_targets())
    print(
        f"Found {found_orthologs} ortholog record(s) across "
        f"{configured_species} target species."
    )
    print(f"Saved {output_path}")
    return 0


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    raise SystemExit(main())
