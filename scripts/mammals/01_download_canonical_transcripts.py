#!/usr/bin/env python3
"""Download canonical protein-coding transcript annotations from Ensembl GTFs."""

import argparse
import csv
import gzip
import logging
import os
import re
import sys
import tempfile
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY_ROOT = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from constants import SPECIES_DICT, SPECIES_ORDER  # noqa: E402

REST_URL = "https://rest.ensembl.org"
FTP_URL = "https://ftp.ensembl.org/pub"
OUTPUT_HEADER = [
    "ensembl_gene_id",
    "gene_name",
    "gene_biotype",
    "ensembl_transcript_id",
    "transcript_biotype",
    "ensembl_protein_id",
    "cds_length",
]
LOGGER = logging.getLogger("canonical_transcripts")


def make_session():
    """Build an HTTP session with retries for transient service errors."""
    retry = Retry(
        total=3,
        backoff_factor=1,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET",),
    )
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retry))
    session.headers.update({"User-Agent": "start-ss mammalian annotations/1.0"})
    session.headers.update({"Accept": "application/json"})
    return session


def get_json(session, url):
    response = session.get(url, timeout=(20, 180))
    response.raise_for_status()
    try:
        return response.json()
    except requests.exceptions.JSONDecodeError as error:
        raise RuntimeError(f"Expected JSON from {url}, received invalid JSON") from error


def normalize_name(value):
    """Normalize a scientific or common name for matching."""
    return re.sub(r"\s+", " ", str(value or "").replace("_", " ").strip()).casefold()


def scientific_name_matches(expected, observed):
    """Match binomials while allowing project metadata to include a subspecies."""
    expected_parts = normalize_name(expected).split()
    observed_parts = normalize_name(observed).split()
    return len(expected_parts) >= 2 and observed_parts[:2] == expected_parts[:2]


def common_name_matches(expected, observed):
    """Match a common name exactly or as a whole phrase within a qualified name."""
    expected_normalized = normalize_name(expected)
    observed_normalized = normalize_name(observed)
    if not expected_normalized or not observed_normalized:
        return False
    return (
        expected_normalized == observed_normalized
        or f" {expected_normalized} " in f" {observed_normalized} "
    )


def find_ensembl_species(species_records, species_id, species_info):
    """Find exactly one REST species record using both configured names."""
    matches = []
    scientific_candidates = []
    for record in species_records:
        display_name = record.get("display_name", "")
        common_name = record.get("common_name", "")
        if not scientific_name_matches(species_info["species"], display_name):
            continue
        scientific_candidates.append(record)
        if common_name_matches(species_info["title"], common_name):
            matches.append(record)

    if len(matches) == 1:
        record = matches[0]
        if not record.get("name") or not record.get("assembly"):
            raise LookupError(
                f"Ensembl species record for {species_id} lacks its production name or assembly"
            )
        return record
    if len(matches) > 1:
        labels = ", ".join(
            f"{item.get('name')} ({item.get('display_name')}; {item.get('common_name')})"
            for item in matches
        )
        raise LookupError(f"Ambiguous Ensembl species match for {species_id}: {labels}")

    if scientific_candidates:
        labels = ", ".join(
            f"{item.get('name')} ({item.get('display_name')}; {item.get('common_name')})"
            for item in scientific_candidates
        )
        raise LookupError(
            f"Scientific name matched for {species_id}, but common name "
            f"{species_info['title']!r} did not; Ensembl candidates: {labels}"
        )
    raise LookupError(
        f"No Ensembl species matched {species_id} "
        f"({species_info['title']}; {species_info['species']})"
    )


def current_release(data_info):
    releases = data_info.get("releases", [])
    numeric_releases = [int(release) for release in releases if str(release).isdigit()]
    if not numeric_releases:
        raise LookupError(f"Ensembl /info/data did not report a usable release: {releases!r}")
    return max(numeric_releases)


def gtf_location(species_record, release):
    """Return the local path and matching FTP URL for a release GTF."""
    production_name = species_record["name"]
    scientific_file_name = species_record["display_name"].replace(" ", "_")
    assembly = species_record["assembly"]
    filename = f"{scientific_file_name}.{assembly}.{release}.gtf.gz"
    url = f"{FTP_URL}/release-{release}/gtf/{production_name}/{filename}"
    path = REPOSITORY_ROOT / "data" / "mammals" / filename
    return path, url


def parse_attributes(attribute_text):
    """Parse the key/value attributes from one GTF record."""
    attributes = {}
    for item in attribute_text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        key, separator, value = item.partition(" ")
        if not separator:
            continue
        value = value.strip().strip('"')
        if key == "tag":
            attributes.setdefault(key, []).append(value)
        else:
            attributes[key] = value
    return attributes


def parse_gtf(gtf_path):
    """Stream a GTF and collect canonical protein-coding transcripts and CDS sizes."""
    transcripts = {}
    with gzip.open(gtf_path, "rt", encoding="utf-8") as gtf_file:
        for line_number, line in enumerate(gtf_file, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed GTF record at line {line_number}: expected 9 columns")

            feature = fields[2]
            if feature not in {"transcript", "CDS"}:
                continue
            attributes = parse_attributes(fields[8])
            transcript_id = attributes.get("transcript_id")
            if not transcript_id:
                continue

            if feature == "transcript":
                tags = [
                    tag
                    for value in attributes.get("tag", [])
                    for tag in value.split(",")
                ]
                if (
                    attributes.get("transcript_biotype") == "protein_coding"
                    and attributes.get("gene_biotype") == "protein_coding"
                    and "Ensembl_canonical" in tags
                ):
                    transcripts[transcript_id] = {
                        "ensembl_gene_id": attributes.get("gene_id", ""),
                        "gene_name": attributes.get("gene_name", ""),
                        "gene_biotype": attributes.get("gene_biotype", ""),
                        "ensembl_transcript_id": transcript_id,
                        "transcript_biotype": attributes.get("transcript_biotype", ""),
                        "ensembl_protein_id": attributes.get("protein_id", ""),
                        "cds_length": 0,
                    }
            elif transcript_id in transcripts:
                try:
                    start, end = int(fields[3]), int(fields[4])
                except ValueError as error:
                    raise ValueError(f"Invalid CDS coordinates at GTF line {line_number}") from error
                if end < start:
                    raise ValueError(f"Invalid CDS interval at GTF line {line_number}: {start}-{end}")
                transcript = transcripts[transcript_id]
                transcript["cds_length"] += end - start + 1
                transcript["ensembl_protein_id"] = (
                    transcript["ensembl_protein_id"] or attributes.get("protein_id", "")
                )

    if not transcripts:
        raise LookupError("GTF contains no canonical protein-coding transcript annotations")

    result = []
    for transcript in transcripts.values():
        if transcript["cds_length"] <= 0:
            raise LookupError(
                "Canonical transcript has no CDS features: "
                f"{transcript['ensembl_transcript_id']}"
            )
        result.append([transcript[field] for field in OUTPUT_HEADER])
    return result


def write_species_file(species_id, rows):
    output_dir = REPOSITORY_ROOT / "data" / "mammals"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"01_canonical_transcripts_{species_id}.tsv"
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{output_path.name}.", dir=output_dir)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as output:
            writer = csv.writer(output, delimiter="\t", lineterminator="\n")
            writer.writerow(OUTPUT_HEADER)
            writer.writerows(rows)
        Path(temporary_name).replace(output_path)
    except Exception:
        Path(temporary_name).unlink(missing_ok=True)
        raise
    return output_path, len(rows)


def process_species(species_id, species_info, species_records, release):
    species_record = find_ensembl_species(species_records, species_id, species_info)
    gtf_path, url = gtf_location(species_record, release)
    LOGGER.info(
        "%s: Ensembl species %s (%s; assembly %s), release %d",
        species_id,
        species_record["name"],
        species_record["display_name"],
        species_record["assembly"],
        release,
    )

    if not gtf_path.is_file():
        relative_path = gtf_path.relative_to(REPOSITORY_ROOT)
        raise FileNotFoundError(
            f"Required Ensembl GTF is missing: {relative_path}. Download it first with "
            f"the wget command documented in docs/2026/1109.mammals/README.md. "
            f"FTP source: {url}"
        )

    rows = parse_gtf(gtf_path)

    output_path, count = write_species_file(species_id, rows)
    LOGGER.info("%s: wrote %d canonical transcripts to %s", species_id, count, output_path)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Download canonical protein-coding transcript annotations from Ensembl GTF files."
    )
    parser.add_argument(
        "species_id",
        help=f"Species ID from constants.py, or __ALL__ for all species ({', '.join(SPECIES_ORDER)})",
    )
    return parser.parse_args()


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    args = parse_args()
    if args.species_id == "__ALL__":
        species_ids = SPECIES_ORDER
    elif args.species_id in SPECIES_DICT:
        species_ids = [args.species_id]
    else:
        LOGGER.error("Unknown species ID %r. Use one of: %s, or __ALL__", args.species_id, ", ".join(SPECIES_ORDER))
        return 2

    session = make_session()
    failures = []
    try:
        species_info_response = get_json(session, f"{REST_URL}/info/species")
        if not isinstance(species_info_response, dict) or not isinstance(
            species_info_response.get("species"), list
        ):
            raise ValueError("Ensembl /info/species response has no species list")
        species_records = species_info_response["species"]
        release = current_release(get_json(session, f"{REST_URL}/info/data"))
    except Exception as error:
        LOGGER.exception("Could not retrieve Ensembl species or release metadata: %s", error)
        return 1

    for species_id in species_ids:
        species_info = SPECIES_DICT.get(species_id)
        if species_info is None:
            message = f"Species ID {species_id!r} is missing from SPECIES_DICT"
            LOGGER.error(message)
            failures.append(message)
            continue
        try:
            process_species(species_id, species_info, species_records, release)
        except Exception as error:
            message = f"{species_id} ({species_info['species']}): {error}"
            LOGGER.exception("Failed to process %s", message)
            failures.append(message)

    if failures:
        LOGGER.error(
            "Completed with %d/%d species failures; outputs may be incomplete:\n- %s",
            len(failures),
            len(species_ids),
            "\n- ".join(failures),
        )
        return 1
    LOGGER.info("Successfully downloaded annotations for %d species", len(species_ids))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
