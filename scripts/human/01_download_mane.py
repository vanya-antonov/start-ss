#!/usr/bin/env python3
"""Download the latest NCBI MANE summary and write MANE Select transcripts."""

from __future__ import annotations

import csv
import gzip
import logging
import re
import sys
import tempfile
from html.parser import HTMLParser
from pathlib import Path
from urllib.error import URLError
from urllib.request import Request, urlopen


FTP_ROOT = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/"
DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "human"
OUTPUT_FILE = DATA_DIR / "mane_select_transcripts.tsv"
OUTPUT_COLUMNS = (
    "symbol",
    "NCBI_GeneID",
    "RefSeq_nuc",
    "Ensembl_nuc",
    "MANE_status",
    "chr",
    "chr_start",
    "chr_end",
    "chr_strand",
)
LOG = logging.getLogger("download_mane")


class DirectoryLinks(HTMLParser):
    """Collect links from the simple NCBI FTP directory index."""

    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "a":
            href = dict(attrs).get("href")
            if href:
                self.links.append(href)


def fetch(url: str) -> bytes:
    """Retrieve a URL with a descriptive user agent and finite timeout."""
    request = Request(url, headers={"User-Agent": "start-ss research script/1.0"})
    with urlopen(request, timeout=60) as response:
        return response.read()


def latest_release_url() -> tuple[str, str]:
    """Find the highest numbered release directory in the official index."""
    parser = DirectoryLinks()
    parser.feed(fetch(FTP_ROOT).decode("utf-8", errors="replace"))
    releases = []
    for link in parser.links:
        match = re.fullmatch(r"release_(\d+(?:\.\d+)*)/", link)
        if match:
            version = tuple(int(part) for part in match.group(1).split("."))
            releases.append((version, match.group(1)))
    if not releases:
        raise RuntimeError(f"No release_* directories found at {FTP_ROOT}")
    _, release = max(releases)
    return release, f"{FTP_ROOT}release_{release}/"


def download_summary(release_url: str, release: str) -> Path:
    """Download this release's compressed summary into data/human."""
    filename = f"MANE.GRCh38.v{release}.summary.txt.gz"
    destination = DATA_DIR / filename
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    payload = fetch(release_url + filename)
    if not payload.startswith(b"\x1f\x8b"):
        raise RuntimeError(f"Downloaded summary is not gzip data: {filename}")
    # Write atomically so an interrupted download cannot leave a partial file.
    with tempfile.NamedTemporaryFile(dir=DATA_DIR, delete=False) as temporary:
        temporary.write(payload)
        temp_path = Path(temporary.name)
    temp_path.replace(destination)
    LOG.info("Downloaded %s (%s bytes)", destination, len(payload))
    return destination


def parse_summary(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read the tab-delimited summary and preserve its actual header names."""
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise RuntimeError("MANE summary has no header row")
        headers = [name.lstrip("#") for name in reader.fieldnames]
        # DictReader keys may have a leading # on the first column; normalize it.
        rows = []
        for raw_row in reader:
            rows.append({key.lstrip("#"): (value or "") for key, value in raw_row.items() if key})
    LOG.info("Summary columns: %s", ", ".join(headers))
    return headers, rows


def filter_mane_select(
    headers: list[str], rows: list[dict[str, str]]
) -> tuple[list[str], list[dict[str, str]]]:
    """Select MANE Select records and map the GRCh38 chromosome label."""
    if "MANE_status" not in headers:
        raise RuntimeError("MANE summary is missing required MANE_status column")
    chromosome_source = next(
        (name for name in ("chr", "GRCh38_chr") if name in headers), None
    )
    selected_columns = [name for name in OUTPUT_COLUMNS if name in headers]
    if "chr" in OUTPUT_COLUMNS and "chr" not in headers and chromosome_source:
        selected_columns.insert(OUTPUT_COLUMNS.index("chr"), "chr")
    if not selected_columns:
        raise RuntimeError("None of the requested output columns are in the summary")

    selected = []
    for row in rows:
        if row.get("MANE_status", "").strip() != "MANE Select":
            continue
        record = {}
        for column in selected_columns:
            source = chromosome_source if column == "chr" else column
            record[column] = row.get(source or column, "")
        selected.append(record)
    if not selected:
        raise RuntimeError("No MANE Select rows found in the summary")
    return selected_columns, selected


def write_table(columns: list[str], rows: list[dict[str, str]], path: Path) -> None:
    """Write a deterministic tab-separated output with a header."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def report_validation(columns: list[str], rows: list[dict[str, str]], release: str) -> None:
    """Log basic data-quality checks and print a few example records."""
    LOG.info("MANE release: %s", release)
    LOG.info("Total MANE Select transcripts: %d", len(rows))
    for column in ("RefSeq_nuc", "NCBI_GeneID"):
        if column not in columns:
            LOG.warning("Cannot check uniqueness: %s was absent", column)
            continue
        values = [row[column] for row in rows if row[column]]
        duplicate_count = len(values) - len(set(values))
        LOG.info("%s unique: %s (duplicates: %d)", column, duplicate_count == 0, duplicate_count)
    for column in columns:
        missing = sum(not row.get(column, "").strip() for row in rows)
        LOG.info("Missing %s: %d/%d", column, missing, len(rows))
    print("Example records:")
    writer = csv.DictWriter(sys.stdout, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows[:5])


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        release, release_url = latest_release_url()
        summary_path = download_summary(release_url, release)
        headers, source_rows = parse_summary(summary_path)
        columns, rows = filter_mane_select(headers, source_rows)
        write_table(columns, rows, OUTPUT_FILE)
        LOG.info("Wrote %s", OUTPUT_FILE)
        report_validation(columns, rows, release)
    except (OSError, URLError, RuntimeError, gzip.BadGzipFile, csv.Error) as error:
        LOG.error("MANE download or processing failed: %s", error)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
