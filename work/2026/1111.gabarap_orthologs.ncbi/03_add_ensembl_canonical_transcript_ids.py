#!/usr/bin/env python3
"""Resolve ortholog symbols and canonical transcripts through Ensembl REST."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen


EXPERIMENT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = EXPERIMENT_DIR.parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "mammals"))
from constants import SPECIES_DICT  # noqa: E402


INPUT_COLUMNS = {"species", "ortholog_gene_symbol", "ensembl_gene_id"}
OUTPUT_COLUMNS = (
    "ensembl_gene_id_resolved",
    "ensembl_canonical_transcript_id",
    "ensembl_lookup_method",
    "ensembl_api_error",
)
ENSEMBL_REST_ROOT = "https://rest.ensembl.org"
ENDPOINT_TEMPLATE = "/lookup/symbol/{species}/{symbol}?expand=1"
CACHE_DIR = EXPERIMENT_DIR / "data" / "ensembl_symbol_lookup_cache"
REQUEST_INTERVAL_SECONDS = 0.1
MAX_ATTEMPTS = 5
TIMEOUT_SECONDS = 60
LOG = logging.getLogger("add_ensembl_canonical_transcript_ids")


class EnsemblRequestError(RuntimeError):
    """An HTTP, connection, timeout, or response decoding failure."""

    def __init__(self, message: str, status_code: int | None = None) -> None:
        super().__init__(message)
        self.status_code = status_code


class SymbolValidationError(ValueError):
    """The symbol lookup response does not match the requested gene symbol."""


@dataclass(frozen=True)
class LookupResult:
    """Normalized outcome for one species and ortholog symbol lookup."""

    resolved_gene_id: str = "-"
    canonical_transcript_id: str = "-"
    lookup_method: str = "-"
    api_error: str = "-"
    category: str = "other_error"
    from_cache: bool = False


def species_name_map() -> dict[str, str]:
    """Map configured scientific names to their Ensembl species names."""
    result: dict[str, str] = {}
    for species_key, details in SPECIES_DICT.items():
        species = str(details.get("species", "")).strip()
        if "ensembl_species_name" not in details:
            raise ValueError(
                f"Species key {species_key!r} in SPECIES_DICT is missing "
                "ensembl_species_name."
            )
        ensembl_name = str(details.get("ensembl_species_name", "")).strip()
        if species and not ensembl_name:
            raise ValueError(
                f"Species {species!r} has an empty ensembl_species_name in "
                "scripts/mammals/constants.py."
            )
        if species:
            result[species.casefold()] = ensembl_name
    return result


def ensembl_species_for(species: str, mapping: dict[str, str]) -> str:
    """Resolve a TSV species label or raise a clear configuration error."""
    key = species.strip().casefold()
    if not key or key not in mapping:
        raise ValueError(
            f"Species {species!r} is not present in SPECIES_DICT; add its "
            "species and ensembl_species_name mapping in scripts/mammals/constants.py."
        )
    return mapping[key]


class EnsemblClient:
    """Rate-limited symbol lookup client with species-symbol keyed JSON cache."""

    def __init__(self, cache_dir: Path) -> None:
        self.cache_dir = cache_dir
        self.next_request_at = 0.0

    def symbol_lookup(
        self, ensembl_species_name: str, gene_symbol: str
    ) -> tuple[dict[str, Any], bool]:
        """Fetch or load the raw symbol lookup response for one key."""
        key = f"{ensembl_species_name}\0{gene_symbol}".encode("utf-8")
        cache_name = hashlib.sha256(key).hexdigest()
        cache_path = self.cache_dir / f"{cache_name}.json"
        if cache_path.is_file():
            try:
                with cache_path.open(encoding="utf-8") as cache_file:
                    payload = json.load(cache_file)
            except (OSError, json.JSONDecodeError) as error:
                raise EnsemblRequestError(
                    f"Invalid cached JSON response: {error}"
                ) from error
            if not isinstance(payload, dict):
                raise EnsemblRequestError("Invalid cached JSON response: expected an object")
            return payload, True

        endpoint = ENDPOINT_TEMPLATE.format(
            species=quote(ensembl_species_name, safe=""),
            symbol=quote(gene_symbol, safe=""),
        )
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
                    "User-Agent": "start-ss Ensembl symbol lookup workflow",
                },
            )
            try:
                with urlopen(request, timeout=TIMEOUT_SECONDS) as response:
                    raw_bytes = response.read()
                try:
                    raw_response = raw_bytes.decode("utf-8")
                except UnicodeDecodeError as error:
                    raise EnsemblRequestError(
                        f"Invalid JSON response encoding: {error}"
                    ) from error
                try:
                    payload = json.loads(raw_response)
                except json.JSONDecodeError as error:
                    raise EnsemblRequestError(
                        f"Invalid JSON response: {error.msg}"
                    ) from error
                if not isinstance(payload, dict):
                    raise EnsemblRequestError(
                        "Invalid JSON response: expected an object"
                    )
                self.cache_dir.mkdir(parents=True, exist_ok=True)
                with cache_path.open("w", encoding="utf-8") as cache_file:
                    json.dump(payload, cache_file, indent=2, ensure_ascii=False)
                    cache_file.write("\n")
                return payload, False
            except HTTPError as error:
                retry_after = error.headers.get("Retry-After")
                description = describe_http_error(error)
                message = f"HTTP {error.code}: {description}"
                retryable = error.code == 429 or 500 <= error.code < 600
                status_code = error.code
            except (URLError, TimeoutError, OSError) as error:
                reason = getattr(error, "reason", error)
                if isinstance(reason, TimeoutError):
                    message = f"Timeout: {reason}"
                else:
                    message = f"Connection error: {reason}"
                retry_after = None
                retryable = True
                status_code = None
            except EnsemblRequestError:
                raise

            if not retryable or attempt == MAX_ATTEMPTS:
                raise EnsemblRequestError(message, status_code=status_code)
            try:
                wait_seconds = float(retry_after) if retry_after else 0.0
            except ValueError:
                wait_seconds = 0.0
            if wait_seconds <= 0:
                wait_seconds = min(2 ** (attempt - 1), 30)
            LOG.warning(
                "Lookup for %s/%s failed (%s); retrying in %.1f seconds",
                ensembl_species_name,
                gene_symbol,
                message,
                wait_seconds,
            )
            time.sleep(wait_seconds)
        raise AssertionError("Unreachable retry loop")


def describe_http_error(error: HTTPError) -> str:
    """Return a concise description from the HTTP response body, if present."""
    try:
        body = error.read().decode("utf-8", errors="replace").strip()
    except OSError:
        body = ""
    if body:
        try:
            payload = json.loads(body)
        except json.JSONDecodeError:
            description = body
        else:
            if isinstance(payload, dict):
                description = str(payload.get("error") or payload.get("message") or body)
            else:
                description = body
        description = " ".join(description.split())
        if len(description) > 240:
            description = description[:237] + "..."
        return description
    return str(error.reason)


def validate_symbol_response(
    payload: dict[str, Any], expected_symbol: str
) -> tuple[str, str]:
    """Validate the resolved symbol and return the gene and transcript IDs."""
    display_name = payload.get("display_name")
    if not isinstance(display_name, str) or not display_name.strip():
        raise SymbolValidationError("Validation error: response has no display_name")
    if display_name.strip().casefold() != expected_symbol.strip().casefold():
        raise SymbolValidationError(
            f"Validation error: expected symbol {expected_symbol!r}, "
            f"but Ensembl returned {display_name!r}"
        )

    resolved_gene_id = payload.get("id")
    if not isinstance(resolved_gene_id, str) or not resolved_gene_id:
        raise EnsemblRequestError("Invalid JSON response: missing top-level gene id")
    canonical_transcript = payload.get("canonical_transcript")
    if canonical_transcript is None:
        transcript_id = "-"
    elif isinstance(canonical_transcript, str):
        transcript_id = canonical_transcript or "-"
    else:
        raise EnsemblRequestError(
            "Invalid JSON response: canonical_transcript was not a string"
        )
    return resolved_gene_id, transcript_id


def read_input_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read the Step 2 TSV and validate columns used by this workflow."""
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
            raise ValueError(
                "Input TSV is missing required column(s): "
                + ", ".join(sorted(missing))
            )
        rows = []
        for line_number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(
                    f"Malformed TSV row at line {line_number}: "
                    "column count differs from header"
                )
            rows.append(row)
        return list(reader.fieldnames), rows


def automatic_paths(gene_symbol: str) -> tuple[Path, Path]:
    """Build Step 2 input and Step 3 output paths from a human gene symbol."""
    safe_symbol = re.sub(
        r"[^A-Za-z0-9._-]+", "_", gene_symbol.strip()
    ).strip("._-")
    if not safe_symbol:
        raise ValueError("Human gene symbol must contain filename-safe characters")
    input_path = EXPERIMENT_DIR / f"02_{safe_symbol}_ncbi_orthologs.tsv"
    output_path = EXPERIMENT_DIR / f"03_{safe_symbol}_ncbi_orthologs.tsv"
    return input_path, output_path


def lookup_result(
    client: EnsemblClient, ensembl_species_name: str, gene_symbol: str
) -> LookupResult:
    """Perform one symbol lookup and normalize its outcome for the table."""
    try:
        payload, from_cache = client.symbol_lookup(ensembl_species_name, gene_symbol)
    except EnsemblRequestError as error:
        category = "not_found" if error.status_code == 404 else "other_error"
        return LookupResult(api_error=str(error), category=category)

    try:
        resolved_gene_id, transcript_id = validate_symbol_response(payload, gene_symbol)
    except SymbolValidationError as error:
        return LookupResult(api_error=str(error), category="validation_error", from_cache=from_cache)
    except EnsemblRequestError as error:
        return LookupResult(api_error=str(error), category="other_error", from_cache=from_cache)

    return LookupResult(
        resolved_gene_id=resolved_gene_id,
        canonical_transcript_id=transcript_id,
        lookup_method="symbol",
        category="success",
        from_cache=from_cache,
    )


def write_output_table(
    output_path: Path, input_columns: list[str], rows: list[dict[str, str]]
) -> None:
    """Write all rows, preserving input columns and appending result columns."""
    columns = list(input_columns)
    for column in OUTPUT_COLUMNS:
        if column not in columns:
            columns.append(column)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(
            output_file, fieldnames=columns, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(
            {column: row.get(column) or "-" for column in columns}
            for row in rows
        )


def process_gene_table(gene_symbol: str, input_path: Path, output_path: Path) -> int:
    """Look up target ortholog symbols and write gene/transcript annotations."""
    print(f"Starting Ensembl symbol lookup for {gene_symbol}.", flush=True)
    print(f"Human gene symbol: {gene_symbol}", flush=True)
    print(f"Input: {input_path}", flush=True)
    print(f"Output: {output_path}", flush=True)
    print(f"Ensembl endpoint: {ENSEMBL_REST_ROOT}{ENDPOINT_TEMPLATE}", flush=True)
    print(f"Cache directory: {CACHE_DIR}", flush=True)

    input_columns, rows = read_input_table(input_path)
    species_mapping = species_name_map()
    row_species: list[tuple[str, str]] = []
    for row in rows:
        species = (row.get("species") or "").strip()
        row_species.append((species, ensembl_species_for(species, species_mapping)))

    client = EnsemblClient(CACHE_DIR)
    results: dict[tuple[str, str, str], LookupResult] = {}
    for row, (species, ensembl_species_name) in zip(rows, row_species):
        ortholog_symbol = (row.get("ortholog_gene_symbol") or "").strip()
        if not ortholog_symbol or ortholog_symbol == "-":
            continue
        key = (species, ensembl_species_name, ortholog_symbol)
        if key in results:
            continue
        if not ensembl_species_name or ensembl_species_name == "-":
            results[key] = LookupResult(
                api_error=(
                    f"Ensembl species name is unavailable for configured species {species!r}"
                ),
                category="other_error",
            )
            continue
        results[key] = lookup_result(client, ensembl_species_name, ortholog_symbol)

    counts = {
        "success": 0,
        "not_found": 0,
        "missing_transcript": 0,
        "other_error": 0,
        "validation_error": 0,
    }
    cache_count = 0
    for (species_label, _, ortholog_symbol), result in results.items():
        if result.category == "success":
            counts["success"] += 1
            if result.canonical_transcript_id == "-":
                counts["missing_transcript"] += 1
        else:
            counts[result.category] += 1
        cache_count += int(result.from_cache)
        print(
            f"{species_label} | {ortholog_symbol}: "
            f"resolved_gene_id={result.resolved_gene_id}; "
            f"canonical_transcript_id={result.canonical_transcript_id}; "
            f"status={result.category}; error={result.api_error}",
            flush=True,
        )

    for row, (species, ensembl_species_name) in zip(rows, row_species):
        ortholog_symbol = (row.get("ortholog_gene_symbol") or "").strip()
        if not ortholog_symbol or ortholog_symbol == "-":
            result = LookupResult(category="skipped")
        else:
            result = results[(species, ensembl_species_name, ortholog_symbol)]
        row["ensembl_gene_id_resolved"] = result.resolved_gene_id
        row["ensembl_canonical_transcript_id"] = result.canonical_transcript_id
        row["ensembl_lookup_method"] = result.lookup_method
        row["ensembl_api_error"] = result.api_error

    write_output_table(output_path, input_columns, rows)
    print("\nAnalysis completed.", flush=True)
    print(f"Successful symbol lookups: {counts['success']}", flush=True)
    print(f"Genes not found (HTTP 404): {counts['not_found']}", flush=True)
    print(f"Missing canonical transcripts: {counts['missing_transcript']}", flush=True)
    print(f"Other API/configuration errors: {counts['other_error']}", flush=True)
    print(f"Symbol validation errors: {counts['validation_error']}", flush=True)
    print(f"Cached responses used: {cache_count}", flush=True)
    print(f"Input rows: {len(rows)}", flush=True)
    print(f"Output: {output_path}", flush=True)
    return 1 if counts["other_error"] or counts["validation_error"] else 0


def main() -> int:
    """Parse the human gene symbol and run the symbol lookup workflow."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "gene_symbol", help="Human gene symbol used in the input and output filenames"
    )
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
