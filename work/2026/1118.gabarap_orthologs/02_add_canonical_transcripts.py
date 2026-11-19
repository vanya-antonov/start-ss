#!/usr/bin/env python3
"""Add Ensembl gene and canonical transcript IDs from Ensembl GTF files."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path


EXPERIMENT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = EXPERIMENT_DIR.parents[2]
DEFAULT_GTF_DIR = PROJECT_ROOT / "data" / "large_files" / "ensembl" / "gtf"
OUTPUT_COLUMNS = ("ensembl_gene", "ensembl_canonical_transcript")

sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "mammals"))
from constants import SPECIES_DICT  # noqa: E402


class GTFFormatError(ValueError):
    """Raised when a GTF record cannot be interpreted safely."""


def parse_attributes(attribute_text: str, line_number: int) -> dict[str, list[str]]:
    """Parse GTF attributes, preserving repeated keys such as ``tag``."""
    attributes: dict[str, list[str]] = defaultdict(list)
    position = 0

    while position < len(attribute_text):
        while position < len(attribute_text) and attribute_text[position] in " ;\t\r\n":
            position += 1
        if position >= len(attribute_text):
            break

        key_start = position
        while position < len(attribute_text) and attribute_text[position] not in " ;\t\r\n":
            position += 1
        key = attribute_text[key_start:position]
        while position < len(attribute_text) and attribute_text[position].isspace():
            position += 1
        if not key or position >= len(attribute_text):
            raise GTFFormatError(f"Malformed attribute at GTF line {line_number}")

        if attribute_text[position] == '"':
            position += 1
            value_chars = []
            while position < len(attribute_text):
                char = attribute_text[position]
                if char == "\\" and position + 1 < len(attribute_text):
                    position += 1
                    value_chars.append(attribute_text[position])
                elif char == '"':
                    position += 1
                    break
                else:
                    value_chars.append(char)
                position += 1
            else:
                raise GTFFormatError(
                    f"Unclosed quoted attribute at GTF line {line_number}"
                )
            value = "".join(value_chars)
        else:
            value_start = position
            while position < len(attribute_text) and attribute_text[position] != ";":
                position += 1
            value = attribute_text[value_start:position].strip()

        while position < len(attribute_text) and attribute_text[position].isspace():
            position += 1
        if position < len(attribute_text):
            if attribute_text[position] != ";":
                raise GTFFormatError(
                    f"Malformed attribute separator at GTF line {line_number}"
                )
            position += 1
        attributes[key].append(value)

    return dict(attributes)


def first_attribute(attributes: dict[str, list[str]], key: str) -> str:
    """Return the first value for an attribute, or an empty string."""
    values = attributes.get(key, [])
    return values[0] if values else ""


def versioned_id(
    attributes: dict[str, list[str]], id_key: str, version_key: str
) -> str:
    """Return an Ensembl ID with its version, whether embedded or separate."""
    identifier = first_attribute(attributes, id_key)
    version = first_attribute(attributes, version_key)
    if identifier and version and not identifier.rsplit(".", maxsplit=1)[-1].isdigit():
        identifier = f"{identifier}.{version}"
    return identifier


def gtf_path_for_species(species_name: str, gtf_dir: Path) -> Path:
    """Find the release 116 GTF using the scientific-name filename prefix."""
    if species_name not in SPECIES_DICT:
        raise ValueError(
            f"Species key {species_name!r} is not present in SPECIES_DICT."
        )
    scientific_name = str(SPECIES_DICT[species_name].get("species", "")).strip()
    if not scientific_name:
        raise ValueError(f"Species key {species_name!r} has no scientific name.")
    filename_prefix = scientific_name.replace(" ", "_")
    filename_prefix = filename_prefix[0].upper() + filename_prefix[1:]
    matches = sorted(gtf_dir.glob(f"{filename_prefix}.*.116.gtf.gz"))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected one release 116 GTF for {species_name!r} matching "
            f"{filename_prefix}.*.116.gtf.gz in {gtf_dir}; found {len(matches)}."
        )
    return matches[0]


def read_input_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read and validate the ortholog TSV."""
    if not path.is_file():
        raise FileNotFoundError(f"Input TSV not found: {path}")
    with path.open(encoding="utf-8-sig", newline="") as input_file:
        reader = csv.DictReader(input_file, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Input TSV is empty or has no header")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("Input TSV contains duplicate column names")
        missing = {"species_name", "ortholog_gene_symbol"} - set(reader.fieldnames)
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


def read_species_gtf(
    gtf_path: Path, requested_symbols: set[str]
) -> tuple[
    dict[str, dict[str, str]],
    dict[str, set[str]],
    dict[str, dict[str, tuple[str, str]]],
]:
    """Index requested gene names and canonical transcript IDs in one GTF pass."""
    genes_by_symbol: dict[str, dict[str, str]] = defaultdict(dict)
    canonical_by_gene: dict[str, set[str]] = defaultdict(set)
    transcript_biotypes: dict[str, dict[str, tuple[str, str]]] = defaultdict(dict)

    try:
        with gzip.open(gtf_path, "rt", encoding="utf-8") as gtf_file:
            for line_number, line in enumerate(gtf_file, start=1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 9:
                    raise GTFFormatError(
                        f"Expected 9 tab-separated columns at GTF line {line_number}; "
                        f"found {len(fields)}"
                    )

                feature = fields[2]
                if feature not in {"gene", "transcript"}:
                    continue
                attributes = parse_attributes(fields[8], line_number)
                gene_name = first_attribute(attributes, "gene_name")
                if gene_name not in requested_symbols:
                    continue
                gene_id = versioned_id(attributes, "gene_id", "gene_version")
                if not gene_id:
                    continue

                if feature == "gene":
                    genes_by_symbol[gene_name][gene_id] = first_attribute(
                        attributes, "gene_biotype"
                    )
                    continue

                tags = {
                    tag
                    for value in attributes.get("tag", [])
                    for tag in value.split(",")
                }
                if "Ensembl_canonical" not in tags:
                    continue
                transcript_id = versioned_id(
                    attributes, "transcript_id", "transcript_version"
                )
                if transcript_id:
                    canonical_by_gene[gene_id].add(transcript_id)
                    transcript_biotypes[gene_id][transcript_id] = (
                        first_attribute(attributes, "gene_biotype"),
                        first_attribute(attributes, "transcript_biotype"),
                    )
    except (OSError, EOFError, UnicodeError) as error:
        raise GTFFormatError(f"Could not read gzipped GTF {gtf_path}: {error}") from error

    return genes_by_symbol, canonical_by_gene, transcript_biotypes


def write_output(
    path: Path, input_columns: list[str], rows: list[dict[str, str]]
) -> None:
    """Write all input rows and append the Ensembl result columns."""
    columns = list(input_columns)
    for column in OUTPUT_COLUMNS:
        if column not in columns:
            columns.append(column)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=columns,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def automatic_paths(gene_name: str) -> tuple[Path, Path]:
    """Build input and output paths from the required gene name."""
    safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", gene_name.strip()).strip("._-")
    if not safe_name:
        raise ValueError("Gene name must contain filename-safe characters")
    return (
        EXPERIMENT_DIR / f"01_{safe_name}_ncbi_orthologs.tsv",
        EXPERIMENT_DIR / f"02_{safe_name}_ncbi_orthologs.tsv",
    )


def process_table(input_path: Path, output_path: Path, gtf_dir: Path) -> int:
    """Resolve gene and canonical transcript IDs for every input row."""
    input_columns, rows = read_input_table(input_path)
    rows_by_species: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        species_name = row["species_name"]
        if species_name not in SPECIES_DICT:
            raise ValueError(
                f"Species key {species_name!r} is not present in SPECIES_DICT."
            )
        rows_by_species[species_name].append(row)

    results: dict[tuple[str, str], tuple[str, str, str]] = {}
    for species_name, species_rows in rows_by_species.items():
        gtf_path = gtf_path_for_species(species_name, gtf_dir)
        symbols = {
            (row.get("ortholog_gene_symbol") or "")
            for row in species_rows
            if (row.get("ortholog_gene_symbol") or "") not in {"", "-"}
        }
        print(f"Scanning {gtf_path} for {len(symbols)} gene symbol(s).", flush=True)
        genes_by_symbol, canonical_by_gene, transcript_biotypes = read_species_gtf(
            gtf_path, symbols
        )

        for row in species_rows:
            symbol = row.get("ortholog_gene_symbol") or ""
            matches = genes_by_symbol.get(symbol, {})
            gene_ids = set(matches)
            if not symbol or symbol == "-" or not gene_ids:
                results[(species_name, symbol)] = ("", "", "missing")
            elif len(gene_ids) > 1:
                results[(species_name, symbol)] = ("", "", "ambiguous")
            else:
                gene_id = next(iter(gene_ids))
                candidate_transcripts = {
                    transcript_id
                    for transcript_id in canonical_by_gene.get(gene_id, set())
                    if transcript_biotypes.get(gene_id, {}).get(transcript_id, ("", ""))[1]
                    == "protein_coding"
                    and transcript_biotypes.get(gene_id, {}).get(transcript_id, ("", ""))[0]
                    in {"", "protein_coding"}
                }
                gene_biotype = matches[gene_id]
                if gene_biotype != "protein_coding":
                    candidate_transcripts = set()
                if len(candidate_transcripts) == 1:
                    results[(species_name, symbol)] = (
                        gene_id,
                        next(iter(candidate_transcripts)),
                        "success",
                    )
                elif len(candidate_transcripts) > 1:
                    results[(species_name, symbol)] = (gene_id, "", "ambiguous_transcript")
                else:
                    results[(species_name, symbol)] = (gene_id, "", "missing_transcript")

    counts = {
        "success": 0,
        "missing": 0,
        "ambiguous": 0,
        "missing_transcript": 0,
        "ambiguous_transcript": 0,
    }
    for row_number, row in enumerate(rows, start=1):
        species_name = row["species_name"]
        symbol = row.get("ortholog_gene_symbol") or ""
        gene_id, transcript_id, status = results[(species_name, symbol)]
        row["ensembl_gene"] = gene_id
        row["ensembl_canonical_transcript"] = transcript_id
        counts[status] += 1
        if status != "success":
            print(
                f"Row {row_number} ({species_name}, {symbol!r}): {status.replace('_', ' ')}.",
                flush=True,
            )

    write_output(output_path, input_columns, rows)
    print("\nAnalysis completed.", flush=True)
    print(f"Successful matches: {counts['success']}", flush=True)
    print(f"Missing genes: {counts['missing']}", flush=True)
    print(f"Ambiguous gene matches: {counts['ambiguous']}", flush=True)
    print(f"Missing canonical transcripts: {counts['missing_transcript']}", flush=True)
    print(f"Ambiguous canonical transcripts: {counts['ambiguous_transcript']}", flush=True)
    print(f"Input rows: {len(rows)}", flush=True)
    print(f"Output: {output_path}", flush=True)
    return 0


def main() -> int:
    """Parse paths and add Ensembl IDs to the GABARAP ortholog table."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "gene_name",
        help="Gene name used to select the experiment input/output filenames",
    )
    parser.add_argument("--input", type=Path, help="Override the input TSV path")
    parser.add_argument("--output", type=Path, help="Override the output TSV path")
    parser.add_argument(
        "--gtf-dir",
        type=Path,
        default=DEFAULT_GTF_DIR,
        help="Directory with Ensembl 116 GTF files",
    )
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        return 2
    args = parser.parse_args()
    try:
        default_input, default_output = automatic_paths(args.gene_name)
        input_path = args.input or default_input
        output_path = args.output or default_output
        return process_table(input_path, output_path, args.gtf_dir)
    except (OSError, UnicodeError, csv.Error, GTFFormatError, ValueError) as error:
        parser.error(f"Could not add canonical transcript IDs: {error}")


if __name__ == "__main__":
    raise SystemExit(main())
