#!/usr/bin/env python3
"""Append Ensembl's canonical transcript ID for each Ensembl gene ID."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import time
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen


EXPERIMENT_DIR = Path(__file__).resolve().parent
INPUT_COLUMNS = {
    "human_gene_symbol",
    "species",
    "ensembl_gene_id",
}
OUTPUT_COLUMN = "ensembl_canonical_transcript_id"
ENSEMBL_REST_ROOT = "https://rest.ensembl.org"
ENDPOINT_TEMPLATE = "/lookup/id/{id}?expand=0"
CACHE_DIR = EXPERIMENT_DIR / "data" / "ensembl_lookup_cache"
REQUEST_INTERVAL_SECONDS = 0.1
MAX_ATTEMPTS = 5
TIMEOUT_SECONDS = 60
LOG = logging.getLogger("add_ensembl_canonical_transcript_ids")


class EnsemblRequestError(RuntimeError):
    """Raised when Ensembl REST lookup fails after retry attempts."""


class EnsemblClient:
    """Rate-limited Ensembl REST client with one raw JSON cache per gene ID."""

    def __init__(self, cache_dir: Path) -> None:
        self.cache_dir = cache_dir
        self.next_request_at = 0.0

    def gene_metadata(self, ensembl_gene_id: str) -> tuple[dict[str, Any], bool]:
        """Return gene lookup metadata and whether it came from cache."""
        cache_path = self.cache_dir / f"{safe_token(ensembl_gene_id)}.json"
        if cache_path.is_file():
            try:
                with cache_path.open(encoding="utf-8") as cache_file:
                    payload = json.load(cache_file)
            except (OSError, json.JSONDecodeError) as error:
                raise EnsemblRequestError(f"Cached response is unreadable: {error}") from error
            return validate_gene_metadata(payload, ensembl_gene_id), True

        endpoint = ENDPOINT_TEMPLATE.format(id=quote(ensembl_gene_id, safe=""))
        url = f"{ENSEMBL_REST_ROOT}{endpoint}"
        for attempt in range(1, MAX_ATTEMPTS + 1):
            delay = self.next_request_at - time.monotonic()
            if delay > 0:
                time.sleep(delay)
            self.next_request_at = time.monotonic() + REQUEST_INTERVAL_SECONDS
            request = Request(
                url,
                headers={
                    "Accept": "application/json",
                    "User-Agent": "start-ss Ensembl canonical transcript workflow",
                },
            )
            try:
                with urlopen(request, timeout=TIMEOUT_SECONDS) as response:
                    payload = json.loads(response.read().decode("utf-8"))
                metadata = validate_gene_metadata(payload, ensembl_gene_id)
                self.cache_dir.mkdir(parents=True, exist_ok=True)
                with cache_path.open("w", encoding="utf-8") as cache_file:
                    json.dump(payload, cache_file, indent=2, ensure_ascii=False)
                    cache_file.write("\n")
                return metadata, False
            except HTTPError as error:
                retry_after = error.headers.get("Retry-After")
                retryable = error.code == 429 or 500 <= error.code < 600
                message = f"HTTP {error.code}"
            except (URLError, TimeoutError, json.JSONDecodeError, ValueError) as error:
                retry_after = None
                retryable = isinstance(error, (URLError, TimeoutError))
                message = str(error)

            if not retryable or attempt == MAX_ATTEMPTS:
                raise EnsemblRequestError(message)
            try:
                wait_seconds = float(retry_after) if retry_after else 0.0
            except ValueError:
                wait_seconds = 0.0
            if wait_seconds <= 0:
                wait_seconds = min(2 ** (attempt - 1), 30)
            LOG.warning(
                "Ensembl lookup for %s failed (%s); retrying in %.1f s",
                ensembl_gene_id,
                message,
                wait_seconds,
            )
            time.sleep(wait_seconds)
        raise AssertionError("Unreachable retry loop")


def safe_token(value: str) -> str:
    """Return a filename-safe token for a stable identifier."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value)


def validate_gene_metadata(payload: Any, requested_id: str) -> dict[str, Any]:
    """Check that an Ensembl lookup response is a gene metadata object."""
    if not isinstance(payload, dict):
        raise ValueError("Ensembl returned JSON that was not an object")
    returned_id = payload.get("id")
    object_type = payload.get("object_type")
    if returned_id and returned_id != requested_id:
        raise ValueError(f"Ensembl returned ID {returned_id!r} for requested ID {requested_id!r}")
    if object_type and object_type != "Gene":
        raise ValueError(f"Ensembl ID {requested_id} refers to {object_type}, not a gene")
    canonical = payload.get("canonical_transcript")
    if canonical is not None and not isinstance(canonical, str):
        raise ValueError("Ensembl canonical_transcript field was not a string")
    return payload


def read_input_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read the Step 2 table and validate columns used by this workflow."""
    if not path.is_file():
        raise FileNotFoundError(f"Input TSV not found: {path}")
    with path.open(encoding="utf-8-sig", newline="") as input_file:
        reader = csv.DictReader(input_file, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Input TSV is empty or has no header")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("Input TSV contains duplicate column names")
        missing = INPUT_COLUMNS - set(reader.fieldnames)
        if missing:
            raise ValueError("Input TSV is missing required column(s): " + ", ".join(sorted(missing)))
        rows = []
        for line_number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Malformed TSV row at line {line_number}: column count differs from header")
            rows.append(row)
        return list(reader.fieldnames), rows


def ensembl_ids_from_rows(rows: list[dict[str, str]]) -> set[str]:
    """Collect distinct Ensembl Gene IDs, allowing semicolon-separated values."""
    return {
        identifier
        for row in rows
        for identifier in split_ensembl_ids(row.get("ensembl_gene_id") or "")
    }


def split_ensembl_ids(value: str) -> list[str]:
    """Split semicolon-separated IDs and discard empty or placeholder values."""
    return [
        identifier
        for raw_identifier in value.split(";")
        if (identifier := raw_identifier.strip()) and identifier != "-"
    ]


def gene_species_map(rows: list[dict[str, str]]) -> dict[str, list[str]]:
    """Associate each Ensembl Gene ID with species labels from input rows."""
    species_by_gene: dict[str, set[str]] = {}
    for row in rows:
        species = (row.get("species") or "").strip() or "Unknown species"
        for gene_id in split_ensembl_ids(row.get("ensembl_gene_id") or ""):
            species_by_gene.setdefault(gene_id, set()).add(species)
    return {gene_id: sorted(species) for gene_id, species in species_by_gene.items()}


def automatic_paths(gene_symbol: str) -> tuple[Path, Path]:
    """Build the Step 2 input and Step 3 output paths from a human gene symbol."""
    safe_symbol = re.sub(r"[^A-Za-z0-9._-]+", "_", gene_symbol.strip()).strip("._-")
    if not safe_symbol:
        raise ValueError("Human gene symbol must contain filename-safe characters")
    input_path = EXPERIMENT_DIR / f"02_{safe_symbol}_ncbi_orthologs.tsv"
    output_path = EXPERIMENT_DIR / f"03_{safe_symbol}_ncbi_orthologs.tsv"
    return input_path, output_path


def write_output_table(
    output_path: Path, input_columns: list[str], rows: list[dict[str, str]]
) -> None:
    """Write all input rows and append or replace the canonical transcript column."""
    columns = list(input_columns)
    if OUTPUT_COLUMN not in columns:
        columns.append(OUTPUT_COLUMN)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(
            {column: row.get(column) or "-" for column in columns}
            for row in rows
        )


def process_gene_table(gene_symbol: str, input_path: Path, output_path: Path) -> int:
    """Retrieve canonical transcript IDs, annotate rows and write the output."""
    print(f"Starting Ensembl canonical transcript lookup for {gene_symbol}.", flush=True)
    print(f"Human gene symbol: {gene_symbol}", flush=True)
    print(f"Input: {input_path}", flush=True)
    print(f"Output: {output_path}", flush=True)
    print(f"Ensembl endpoint: {ENSEMBL_REST_ROOT}{ENDPOINT_TEMPLATE}", flush=True)
    print(f"Cache directory: {CACHE_DIR}", flush=True)
    print("Reading input TSV...", flush=True)
    input_columns, rows = read_input_table(input_path)

    gene_ids = sorted(ensembl_ids_from_rows(rows))
    species_by_gene = gene_species_map(rows)
    missing_id_rows = sum(
        not split_ensembl_ids(row.get("ensembl_gene_id") or "") for row in rows
    )
    print(f"Input rows: {len(rows)}", flush=True)
    print(f"Unique Ensembl Gene IDs: {len(gene_ids)}", flush=True)
    print(f"Rows without an Ensembl Gene ID: {missing_id_rows}", flush=True)

    client = EnsemblClient(CACHE_DIR)
    canonical_by_gene: dict[str, str] = {}
    failures: dict[str, str] = {}
    cache_count = 0
    request_count = 0

    for index, gene_id in enumerate(gene_ids, start=1):
        species = ", ".join(species_by_gene.get(gene_id, [])) or "Unknown species"
        print(f"[{index}/{len(gene_ids)}] Processing {species}: {gene_id}", flush=True)
        cache_file = CACHE_DIR / f"{safe_token(gene_id)}.json"
        if cache_file.is_file():
            print("  Using cached response...", flush=True)
        else:
            print("  Fetching gene metadata from Ensembl REST API...", flush=True)
            request_count += 1
        try:
            metadata, from_cache = client.gene_metadata(gene_id)
            if from_cache:
                cache_count += 1
            canonical_id = metadata.get("canonical_transcript") or ""
            canonical_by_gene[gene_id] = canonical_id
            if canonical_id:
                print(f"  Canonical transcript: {canonical_id}", flush=True)
            else:
                print("  No canonical transcript available.", flush=True)
        except (EnsemblRequestError, ValueError) as error:
            failures[gene_id] = str(error)
            canonical_by_gene[gene_id] = ""
            print(f"  API error or invalid Ensembl Gene ID: {error}", flush=True)

    for row in rows:
        ids = [
            canonical_by_gene.get(gene_id.strip(), "")
            for gene_id in split_ensembl_ids(row.get("ensembl_gene_id") or "")
        ]
        row[OUTPUT_COLUMN] = (
            ";".join(dict.fromkeys(transcript_id for transcript_id in ids if transcript_id))
            or "-"
        )

    print(f"Writing output TSV ({len(rows)} rows)...", flush=True)
    write_output_table(output_path, input_columns, rows)
    print(f"Output written: {output_path}", flush=True)
    print(f"Rows written: {len(rows)}", flush=True)

    no_canonical_count = sum(
        1 for gene_id, transcript_id in canonical_by_gene.items()
        if not transcript_id and gene_id not in failures
    )
    print("\nAnalysis completed.", flush=True)
    print(f"Gene: {gene_symbol}", flush=True)
    print(f"Input rows: {len(rows)}", flush=True)
    print(f"Unique Ensembl Gene IDs: {len(gene_ids)}", flush=True)
    print(f"Canonical transcripts retrieved: {sum(bool(value) for value in canonical_by_gene.values())}", flush=True)
    print(f"Genes without canonical transcript: {no_canonical_count}", flush=True)
    print(f"Genes without Ensembl Gene ID: {missing_id_rows}", flush=True)
    print(f"API failures: {len(failures)}", flush=True)
    print(f"Cached responses used: {cache_count}", flush=True)
    print(f"New API requests: {request_count}", flush=True)
    print(f"Output: {output_path}", flush=True)
    return 1 if failures else 0


def main() -> int:
    """Parse the human gene symbol and run the canonical transcript workflow."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gene_symbol", help="Human gene symbol used in the input and output filenames")
    args = parser.parse_args()
    try:
        input_path, output_path = automatic_paths(args.gene_symbol)
        return process_gene_table(args.gene_symbol, input_path, output_path)
    except (OSError, ValueError, csv.Error) as error:
        LOG.error("Could not retrieve canonical transcript IDs: %s", error)
        return 1


if __name__ == "__main__":
    logging.basicConfig(level=logging.ERROR, format="%(levelname)s: %(message)s")
    raise SystemExit(main())
