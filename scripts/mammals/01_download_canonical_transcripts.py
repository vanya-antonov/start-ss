#!/usr/bin/env python3
"""Download canonical protein-coding transcript annotations from Ensembl BioMart."""

import csv
import argparse
import logging
import os
import sys
import tempfile
from pathlib import Path
import xml.etree.ElementTree as ET

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

SCRIPT_DIR = Path(__file__).resolve().parent
REPOSITORY_ROOT = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from constants import SPECIES_DICT, SPECIES_ORDER  # noqa: E402

BIOMART_URL = "https://www.ensembl.org/biomart/martservice"
MART_NAME = "ENSEMBL_MART_ENSEMBL"
OUTPUT_FIELDS = [
    ("ensembl_gene_id", "Gene stable ID"),
    ("external_gene_name", "Gene name"),
    ("gene_biotype", "Gene biotype"),
    ("ensembl_transcript_id", "Transcript stable ID"),
    ("transcript_biotype", "Transcript type"),
    ("ensembl_peptide_id", "Protein stable ID"),
    ("cds_length", "CDS Length"),
]
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
    """Build an HTTP session with bounded retries for transient API failures."""
    retry = Retry(
        total=3,
        backoff_factor=1,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET", "POST"),
    )
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retry))
    session.headers.update({"User-Agent": "start-ss mammalian annotations/1.0"})
    return session


def get_metadata(session, **params):
    """Fetch a BioMart metadata table and return non-empty tab-separated rows."""
    response = session.get(BIOMART_URL, params=params, timeout=(20, 180))
    response.raise_for_status()
    lines = [line for line in response.text.splitlines() if line.strip()]
    if not lines or "Query ERROR" in response.text:
        raise RuntimeError(f"BioMart returned invalid metadata for {params}: {response.text[:300]}")
    return list(csv.reader(lines, delimiter="\t"))


def find_dataset(dataset_rows, scientific_name):
    """Match a species to a dataset using its scientific name, not its UCSC ID."""
    name_parts = scientific_name.casefold().split()
    if len(name_parts) < 2:
        raise ValueError(f"Expected a binomial scientific name: {scientific_name}")
    # Ensembl dataset labels generally use genus and species, even for a
    # trinomial subspecies name supplied in project metadata.
    required_tokens = name_parts[:2]
    matches = []
    for row in dataset_rows:
        if len(row) < 2:
            continue
        searchable = " ".join(row[1:]).casefold()
        if all(token in searchable for token in required_tokens):
            matches.append((row[0], row[1]))
    if len(matches) != 1:
        details = ", ".join(f"{name} ({label})" for name, label in matches) or "none"
        raise LookupError(
            f"Expected one Ensembl dataset for {scientific_name}; found {len(matches)}: {details}"
        )
    return matches[0]


def metadata_names(rows):
    """Return first-column names from BioMart attribute/filter metadata."""
    return {row[0].strip() for row in rows if row and row[0].strip()}


def choose_filter(available_filters, candidates, description):
    for name in candidates:
        if name in available_filters:
            return name
    raise LookupError(f"BioMart dataset does not expose a {description} filter")


def build_query(dataset, attributes, gene_biotype_filter, canonical_filter):
    query = ET.Element(
        "Query",
        {
            "virtualSchemaName": "default",
            "formatter": "TSV",
            "header": "0",
            "uniqueRows": "1",
            "count": "",
            "datasetConfigVersion": "0.6",
        },
    )
    dataset_element = ET.SubElement(query, "Dataset", {"name": dataset, "interface": "default"})
    ET.SubElement(dataset_element, "Filter", {"name": gene_biotype_filter, "value": "protein_coding"})
    ET.SubElement(dataset_element, "Filter", {"name": canonical_filter, "value": "1"})
    for attribute, _ in attributes:
        ET.SubElement(dataset_element, "Attribute", {"name": attribute})
    return ET.tostring(query, encoding="unicode")


def parse_result(text, requested_attributes, available_attributes):
    """Validate the BioMart response and expand unavailable optional columns."""
    if "Query ERROR" in text or "ERROR:" in text[:500]:
        raise RuntimeError(f"BioMart query failed: {text[:500]}")
    rows = list(csv.reader(text.splitlines(), delimiter="\t"))
    if not rows:
        raise RuntimeError("BioMart returned an empty response")
    expected = [attribute for attribute, _ in requested_attributes]
    if len(rows) == 0:
        raise RuntimeError("BioMart query returned no canonical protein-coding transcripts")

    optional_indices = {
        index for index, (attribute, _) in enumerate(requested_attributes)
        if attribute not in available_attributes
    }
    result = []
    for line_number, row in enumerate(rows, start=1):
        if len(row) != len(expected):
            raise RuntimeError(
                f"Malformed BioMart row {line_number}: expected {len(expected)} columns, got {len(row)}"
            )
        for index in optional_indices:
            row[index] = ""
        if not row[0] or not row[3]:
            raise RuntimeError(f"BioMart row {line_number} is missing a gene or transcript ID")
        result.append(row)
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


def download_species(session, species_id, species_info, dataset_rows):
    scientific_name = species_info["species"]
    dataset, display_name = find_dataset(dataset_rows, scientific_name)
    LOGGER.info("%s: matched %s (%s)", species_id, dataset, display_name)

    attributes_rows = get_metadata(session, type="attributes", dataset=dataset, mart=MART_NAME)
    filters_rows = get_metadata(session, type="filters", dataset=dataset, mart=MART_NAME)
    available_attributes = metadata_names(attributes_rows)
    available_filters = metadata_names(filters_rows)

    required_attributes = [
        ("ensembl_gene_id", "Ensembl gene ID"),
        ("external_gene_name", "gene name"),
        ("gene_biotype", "gene biotype"),
        ("ensembl_transcript_id", "Ensembl transcript ID"),
        ("transcript_biotype", "transcript biotype"),
    ]
    missing = [name for name, _ in required_attributes if name not in available_attributes]
    if missing:
        raise LookupError(f"Dataset {dataset} lacks required attributes: {', '.join(missing)}")

    attributes = [
        (name, label)
        for name, label in OUTPUT_FIELDS
        if name in available_attributes or name in {item[0] for item in required_attributes}
    ]
    gene_biotype_filter = choose_filter(
        available_filters, ("biotype", "gene_biotype"), "protein-coding gene biotype"
    )
    canonical_filter = choose_filter(
        available_filters, ("transcript_is_canonical",), "canonical transcript"
    )
    query = build_query(dataset, attributes, gene_biotype_filter, canonical_filter)
    response = session.get(BIOMART_URL, params={"query": query}, timeout=(20, 600))
    response.raise_for_status()
    parsed = parse_result(response.text, attributes, available_attributes)

    # Place each value in the stable seven-column output schema. Optional
    # fields absent from a dataset remain empty instead of shifting columns.
    rows = []
    attribute_indices = {name: index for index, (name, _) in enumerate(attributes)}
    for source_row in parsed:
        rows.append([
            source_row[attribute_indices[name]] if name in attribute_indices else ""
            for name, _ in OUTPUT_FIELDS
        ])
    path, count = write_species_file(species_id, rows)
    LOGGER.info("%s: wrote %d canonical transcripts to %s", species_id, count, path)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Download canonical protein-coding transcript annotations from Ensembl BioMart."
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
        dataset_rows = get_metadata(
            session, type="datasets", mart=MART_NAME, virtualSchema="default"
        )
    except Exception:
        LOGGER.exception("Could not retrieve the Ensembl BioMart dataset registry")
        return 1

    for species_id in species_ids:
        species_info = SPECIES_DICT.get(species_id)
        if species_info is None:
            message = f"Species ID {species_id!r} is missing from SPECIES_DICT"
            LOGGER.error(message)
            failures.append(message)
            continue
        try:
            download_species(session, species_id, species_info, dataset_rows)
        except Exception as error:
            message = f"{species_id} ({species_info['species']}): {error}"
            LOGGER.exception("Failed to download %s", message)
            failures.append(message)

    if failures:
        LOGGER.error(
            "Completed with %d/%d species failures; outputs may be incomplete:\n- %s",
            len(failures),
            len(species_ids),
            "\n- ".join(failures),
        )
        return 1
    LOGGER.info("Successfully downloaded canonical transcript annotations for %d species", len(species_ids))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
