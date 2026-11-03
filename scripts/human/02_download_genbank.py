#!/usr/bin/env python3
"""Download versioned MANE Select transcript records from NCBI EFetch."""

from __future__ import annotations

import csv
import gzip
import logging
import os
import re
import tempfile
import time
from io import StringIO
from pathlib import Path
from typing import TextIO

import requests
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "human"
INPUT_FILE = DATA_DIR / "mane_select_transcripts.tsv"
OUTPUT_FILE = DATA_DIR / "mane_transcripts.gb.gz"
FAILURES_FILE = DATA_DIR / "genbank_download_failures.tsv"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BATCH_SIZE = 150
MAX_ATTEMPTS = 5
LOG = logging.getLogger("download_genbank")
ACCESSION_VERSION = re.compile(r"^[A-Za-z]{1,6}_?[A-Za-z0-9]+\.\d+$")


def read_accessions(path: Path) -> list[str]:
    """Read unique, versioned RefSeq accessions from the Step 1 TSV."""
    if not path.is_file():
        raise FileNotFoundError(f"Input MANE table not found: {path}")
    accessions: list[str] = []
    seen: set[str] = set()
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "RefSeq_nuc" not in reader.fieldnames:
            raise ValueError("Input TSV must contain a RefSeq_nuc column")
        for line_number, row in enumerate(reader, start=2):
            accession = (row.get("RefSeq_nuc") or "").strip()
            if not accession:
                raise ValueError(f"Missing RefSeq_nuc at input line {line_number}")
            if not ACCESSION_VERSION.fullmatch(accession):
                raise ValueError(
                    f"Expected a versioned RefSeq accession at line {line_number}: {accession!r}"
                )
            if accession in seen:
                LOG.warning("Duplicate input accession ignored: %s", accession)
                continue
            seen.add(accession)
            accessions.append(accession)
    if not accessions:
        raise ValueError("Input TSV contains no RefSeq accessions")
    LOG.info("Read %d unique versioned accessions", len(accessions))
    return accessions


def request_parameters() -> dict[str, str]:
    """Build EFetch parameters, optionally using NCBI credentials from env."""
    params = {"db": "nuccore", "rettype": "gb", "retmode": "text", "tool": "start_ss"}
    email = os.environ.get("NCBI_EMAIL", "").strip()
    api_key = os.environ.get("NCBI_API_KEY", "").strip()
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    return params


def wait_for_rate_limit(next_request_at: float, interval: float) -> float:
    """Wait until a request is allowed and return the next permitted time."""
    now = time.monotonic()
    if now < next_request_at:
        time.sleep(next_request_at - now)
    return time.monotonic() + interval


def fetch_batch(
    session: requests.Session,
    accessions: list[str],
    params: dict[str, str],
    interval: float,
    next_request_at: float,
) -> tuple[list[SeqRecord], float]:
    """POST a batch to EFetch with rate limiting and transient-error retries."""
    request_data = {**params, "id": ",".join(accessions)}
    for attempt in range(1, MAX_ATTEMPTS + 1):
        next_request_at = wait_for_rate_limit(next_request_at, interval)
        try:
            response = session.post(EFETCH_URL, data=request_data, timeout=(20, 180))
            if response.status_code == 429 or response.status_code >= 500:
                response.raise_for_status()
            response.raise_for_status()
            try:
                records = list(SeqIO.parse(StringIO(response.text), "genbank"))
            except Exception as error:
                raise RuntimeError(f"Could not parse EFetch GenBank response: {error}") from error
            if not records:
                raise RuntimeError("EFetch returned no parseable GenBank records")
            return records, next_request_at
        except (requests.RequestException, RuntimeError) as error:
            if attempt == MAX_ATTEMPTS:
                raise RuntimeError(
                    f"Batch failed after {MAX_ATTEMPTS} attempts: {error}"
                ) from error
            delay = min(2 ** (attempt - 1), 30)
            LOG.warning(
                "Batch request attempt %d/%d failed (%s); retrying in %d seconds",
                attempt,
                MAX_ATTEMPTS,
                error,
                delay,
            )
            time.sleep(delay)
    raise AssertionError("Unreachable retry loop")


def accession_of(record: SeqRecord) -> str:
    """Return the record identifier including its accession version."""
    return record.id.strip()


def save_records(
    accessions: list[str], output_handle: TextIO, params: dict[str, str]
) -> tuple[set[str], dict[str, str]]:
    """Download batches and write only accession/version-exact records."""
    requested = set(accessions)
    downloaded: set[str] = set()
    failures: dict[str, str] = {}
    batches = [accessions[i : i + BATCH_SIZE] for i in range(0, len(accessions), BATCH_SIZE)]
    interval = 0.11 if params.get("api_key") else 0.34
    next_request_at = 0.0

    with requests.Session() as session:
        for batch_number, batch in enumerate(batches, start=1):
            try:
                records, next_request_at = fetch_batch(
                    session, batch, params, interval, next_request_at
                )
            except RuntimeError as error:
                for accession in batch:
                    failures[accession] = str(error)
                LOG.error("Batch %d/%d failed: %s", batch_number, len(batches), error)
                continue

            seen_in_batch: set[str] = set()
            for record in records:
                accession = accession_of(record)
                if accession not in batch:
                    LOG.warning("Ignoring unexpected EFetch record: %s", accession)
                    continue
                if accession in seen_in_batch or accession in downloaded:
                    LOG.warning("Ignoring duplicate EFetch record: %s", accession)
                    continue
                seen_in_batch.add(accession)
                SeqIO.write(record, output_handle, "genbank")
                downloaded.add(accession)

            for accession in batch:
                if accession not in seen_in_batch:
                    failures[accession] = (
                        "No exact accession.version record returned by EFetch"
                    )
            LOG.info(
                "Batch %d/%d complete: %d/%d records downloaded; %d failed/missing so far",
                batch_number,
                len(batches),
                len(downloaded),
                len(accessions),
                len(failures),
            )

    # Keep the report strictly about expected inputs, including any not returned.
    for accession in requested - downloaded:
        failures.setdefault(accession, "Accession was not downloaded")
    return downloaded, failures


def validate_genbank(path: Path, expected: set[str]) -> tuple[int, int]:
    """Parse the staged file and check exact versions, uniqueness, and CDS."""
    seen: set[str] = set()
    cds_records = 0
    count = 0
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            accession = accession_of(record)
            if accession not in expected:
                raise RuntimeError(f"Unexpected accession/version in output: {accession}")
            if accession in seen:
                raise RuntimeError(f"Duplicate accession/version in output: {accession}")
            seen.add(accession)
            count += 1
            if any(feature.type == "CDS" for feature in record.features):
                cds_records += 1
    LOG.info("Validated accession.version for all %d records", count)
    LOG.info("Unique GenBank records: %d", len(seen))
    LOG.info("Records with an annotated CDS feature: %d", cds_records)
    return count, cds_records


def write_failures(path: Path, failures: dict[str, str]) -> None:
    """Write a failure report only when failures exist; remove stale reports."""
    if not failures:
        path.unlink(missing_ok=True)
        LOG.info("No failed or missing accessions; failure report not created")
        return

    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as temporary:
        writer = csv.writer(temporary, delimiter="\t", lineterminator="\n")
        writer.writerow(("RefSeq_nuc", "reason"))
        for accession in sorted(failures):
            writer.writerow((accession, failures[accession]))
        temp_path = Path(temporary.name)
    temp_path.replace(path)


def download_all(accessions: list[str]) -> None:
    """Stage, validate, and atomically publish the GenBank and failure files."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    params = request_parameters()
    if not params.get("email"):
        LOG.warning("NCBI_EMAIL is unset; provide an email via the NCBI_EMAIL environment variable")
    if params.get("api_key"):
        LOG.info("Using NCBI API key from NCBI_API_KEY")

    with tempfile.NamedTemporaryFile(dir=DATA_DIR, delete=False) as temporary_file:
        staged_output = Path(temporary_file.name)

    try:
        with gzip.open(staged_output, "wt", encoding="utf-8", newline="") as temporary:
            downloaded, failures = save_records(accessions, temporary, params)
        with staged_output.open("rb") as temporary_file:
            os.fsync(temporary_file.fileno())
    except Exception:
        staged_output.unlink(missing_ok=True)
        raise

    try:
        count, _ = validate_genbank(staged_output, set(accessions))
        if count != len(downloaded):
            raise RuntimeError(
                f"Staged record count ({count}) differs from downloaded count ({len(downloaded)})"
            )
        if set(accessions) - downloaded:
            LOG.warning("Missing or failed accessions: %d", len(set(accessions) - downloaded))
        os.replace(staged_output, OUTPUT_FILE)
    except Exception:
        staged_output.unlink(missing_ok=True)
        raise

    write_failures(FAILURES_FILE, failures)
    LOG.info("Successfully downloaded records: %d/%d", count, len(accessions))
    LOG.info("Wrote GenBank records to %s", OUTPUT_FILE)
    if failures:
        LOG.info("Wrote failure report (%d accessions) to %s", len(failures), FAILURES_FILE)


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
    try:
        accessions = read_accessions(INPUT_FILE)
        download_all(accessions)
    except (OSError, ValueError, RuntimeError, csv.Error) as error:
        LOG.error("GenBank download failed: %s", error)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
