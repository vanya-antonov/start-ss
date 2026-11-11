#!/usr/bin/env python3
"""Append the longest annotated protein-coding transcript for each gene."""

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
DATA_DIR = EXPERIMENT_DIR / "data"
API_ROOT = "https://api.ncbi.nlm.nih.gov/datasets/v2"
REQUEST_INTERVAL_SECONDS = 0.21  # NCBI Datasets default rate limit: 5 requests/s.
MAX_ATTEMPTS = 5
TIMEOUT_SECONDS = 60
LOG = logging.getLogger("add_longest_protein_transcript")
REQUIRED_COLUMNS = {
    "taxon_id",
    "ortholog_gene_id",
}
TRANSCRIPT_COLUMNS = [
    "transcript_accession",
    "transcript_length",
    "protein_accession",
    "protein_length",
    "transcript_selection_method",
]


class NCBIRequestError(RuntimeError):
    """Raised when an NCBI request fails after retry attempts."""


class NCBIClient:
    """Rate-limited NCBI Datasets JSON client with raw-response caching."""

    def __init__(self, data_dir: Path) -> None:
        self.data_dir = data_dir
        self.next_request_at = 0.0

    def gene_product_report(self, taxon_id: str, gene_id: str) -> dict[str, Any]:
        """Retrieve a gene's Datasets product report by NCBI GeneID."""
        cache_path = self.data_dir / f"gene_product_taxon_{safe_token(taxon_id)}_gene_{safe_token(gene_id)}.json"
        if cache_path.exists():
            with cache_path.open(encoding="utf-8") as raw_file:
                return json.load(raw_file)

        endpoint = f"/gene/id/{quote(gene_id, safe='')}/product_report"
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
                    "User-Agent": "start-ss NCBI transcript selection workflow",
                },
            )
            try:
                with urlopen(request, timeout=TIMEOUT_SECONDS) as response:
                    payload = json.loads(response.read().decode("utf-8"))
                self.data_dir.mkdir(parents=True, exist_ok=True)
                with cache_path.open("w", encoding="utf-8") as raw_file:
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
                raise NCBIRequestError(message)
            try:
                wait_seconds = float(retry_after) if retry_after else 0.0
            except ValueError:
                wait_seconds = 0.0
            if wait_seconds <= 0:
                wait_seconds = min(2 ** (attempt - 1), 30)
            LOG.warning(
                "Gene %s (taxon %s), attempt %d/%d failed: %s; retrying in %.1f s",
                gene_id,
                taxon_id,
                attempt,
                MAX_ATTEMPTS,
                message,
                wait_seconds,
            )
            time.sleep(wait_seconds)
        raise AssertionError("Unreachable retry loop")


def safe_token(value: str) -> str:
    """Return a filename-safe token for a taxon or GeneID."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value)


def reports_from_response(payload: dict[str, Any]) -> list[dict[str, Any]]:
    """Return product report gene records from the Datasets response wrapper."""
    reports = payload.get("reports")
    if not isinstance(reports, list):
        raise ValueError("NCBI product report did not contain a reports array")
    return [item.get("gene", item) for item in reports if isinstance(item, dict)]


def parse_positive_length(value: Any) -> int | None:
    """Parse a positive annotated length, rejecting absent and malformed values."""
    try:
        length = int(value)
    except (TypeError, ValueError):
        return None
    return length if length > 0 else None


def select_longest_protein_transcript(
    product_report: dict[str, Any], expected_taxon_id: str, expected_gene_id: str
) -> dict[str, str] | None:
    """Select a protein-coding transcript by protein length and stable tie-breaks."""
    candidates: list[dict[str, Any]] = []
    for gene in reports_from_response(product_report):
        if str(gene.get("geneId", "")) != expected_gene_id:
            continue
        if str(gene.get("taxId", "")) != expected_taxon_id:
            continue
        transcripts = gene.get("transcripts")
        if not isinstance(transcripts, list):
            continue
        for transcript in transcripts:
            if not isinstance(transcript, dict):
                continue
            transcript_type = transcript.get("type", "")
            if transcript_type not in {"PROTEIN_CODING", "PROTEIN_CODING_MODEL"}:
                continue
            protein = transcript.get("protein")
            if not isinstance(protein, dict):
                continue
            transcript_accession = str(transcript.get("accessionVersion") or "")
            protein_accession = str(protein.get("accessionVersion") or "")
            protein_length = parse_positive_length(protein.get("length"))
            transcript_length = parse_positive_length(transcript.get("length"))
            if not transcript_accession or not protein_accession or protein_length is None:
                continue
            candidates.append({
                "transcript_accession": transcript_accession,
                "transcript_length": str(transcript_length) if transcript_length else "",
                "protein_accession": protein_accession,
                "protein_length": protein_length,
            })

    if not candidates:
        return None

    def selection_key(candidate: dict[str, Any]) -> tuple[int, int, str]:
        transcript_accession = candidate["transcript_accession"]
        if transcript_accession.startswith("NM_"):
            accession_priority = 0
        else:
            accession_priority = 1
        return (
            -candidate["protein_length"],
            accession_priority,
            transcript_accession,
        )

    selected = min(candidates, key=selection_key)
    return {
        "transcript_accession": selected["transcript_accession"],
        "transcript_length": selected["transcript_length"],
        "protein_accession": selected["protein_accession"],
        "protein_length": str(selected["protein_length"]),
        "transcript_selection_method": "longest_protein",
    }


def automatic_output_path(input_path: Path) -> Path:
    """Create the requested default output filename beside the input file."""
    if re.match(r"^\d+_", input_path.name):
        output_name = re.sub(r"^\d+_", "02_", input_path.name, count=1)
    else:
        output_name = f"{input_path.stem}_with_transcripts{input_path.suffix}"
    return input_path.with_name(output_name)


def read_input_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a TSV and validate required ortholog identifier columns."""
    with path.open(encoding="utf-8-sig", newline="") as input_file:
        reader = csv.DictReader(input_file, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Input TSV has no header")
        missing = REQUIRED_COLUMNS - set(reader.fieldnames)
        if missing:
            raise ValueError("Input TSV is missing columns: " + ", ".join(sorted(missing)))
        rows = list(reader)
        return list(reader.fieldnames), rows


def write_output_table(
    path: Path, input_columns: list[str], rows: list[dict[str, str]]
) -> None:
    """Write the original columns and appended transcript metadata."""
    columns = input_columns + [name for name in TRANSCRIPT_COLUMNS if name not in input_columns]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def process_table(input_path: Path, output_path: Path) -> int:
    """Annotate unique genes, preserve all rows, and write the output TSV."""
    input_columns, rows = read_input_table(input_path)
    data_dir = EXPERIMENT_DIR / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    client = NCBIClient(data_dir)

    gene_keys = sorted({
        (str(row.get("taxon_id", "")).strip(), str(row.get("ortholog_gene_id", "")).strip())
        for row in rows
        if str(row.get("taxon_id", "")).strip()
        and str(row.get("ortholog_gene_id", "")).strip()
        and str(row.get("ortholog_gene_id", "")).strip() != "-"
    })
    selected_by_gene: dict[tuple[str, str], dict[str, str] | None] = {}
    failures: dict[tuple[str, str], str] = {}

    for index, (taxon_id, gene_id) in enumerate(gene_keys, start=1):
        print(f"[{index}/{len(gene_keys)}] NCBI GeneID {gene_id} (taxon {taxon_id})")
        try:
            response = client.gene_product_report(taxon_id, gene_id)
            selected_by_gene[(taxon_id, gene_id)] = select_longest_protein_transcript(
                response, taxon_id, gene_id
            )
        except (NCBIRequestError, OSError, json.JSONDecodeError, ValueError) as error:
            failures[(taxon_id, gene_id)] = str(error)
            selected_by_gene[(taxon_id, gene_id)] = None
            LOG.error("Failed for GeneID %s (taxon %s): %s", gene_id, taxon_id, error)

    for row in rows:
        key = (str(row.get("taxon_id", "")).strip(), str(row.get("ortholog_gene_id", "")).strip())
        selected = selected_by_gene.get(key)
        if selected:
            row.update(selected)
        else:
            row.update({column: "" for column in TRANSCRIPT_COLUMNS})
    write_output_table(output_path, input_columns, rows)

    selected_count = sum(value is not None for value in selected_by_gene.values())
    no_transcript_count = sum(
        value is None and key not in failures for key, value in selected_by_gene.items()
    )
    print("\nSummary")
    print(f"Input rows: {len(rows)}")
    print(f"Unique genes processed: {len(gene_keys)}")
    print(f"Genes with selected transcript: {selected_count}")
    print(f"Genes without suitable transcript: {no_transcript_count}")
    print(f"API failures: {len(failures)}")
    print(f"Output: {output_path}")
    if no_transcript_count:
        without_transcript = [
            f"taxon {key[0]}, GeneID {key[1]}"
            for key, value in selected_by_gene.items()
            if value is None and key not in failures
        ]
        print("No suitable transcript: " + "; ".join(without_transcript))
    return 1 if failures else 0


def main() -> int:
    """Parse arguments and execute transcript annotation."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_tsv", type=Path, help="Ortholog input table in TSV format")
    parser.add_argument("output_tsv", nargs="?", type=Path, help="Optional output TSV path")
    args = parser.parse_args()

    input_path = args.input_tsv.expanduser()
    output_path = args.output_tsv.expanduser() if args.output_tsv else automatic_output_path(input_path)
    try:
        return process_table(input_path, output_path)
    except (OSError, ValueError, csv.Error) as error:
        LOG.error("Could not process input: %s", error)
        return 1


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    raise SystemExit(main())
