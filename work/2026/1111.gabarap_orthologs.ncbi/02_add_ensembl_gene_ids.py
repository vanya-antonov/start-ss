#!/usr/bin/env python3
"""Map NCBI Gene IDs to Ensembl Gene IDs using NCBI's local mapping file."""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import re
from pathlib import Path


EXPERIMENT_DIR = Path(__file__).resolve().parent
REPOSITORY_ROOT = EXPERIMENT_DIR.parents[2]
MAPPING_FILE = REPOSITORY_ROOT / "data" / "large_files" / "ncbi" / "gene2ensembl.gz"
REQUIRED_INPUT_COLUMNS = {"taxon_id", "ortholog_gene_id"}
MAPPING_COLUMNS = ("tax_id", "GeneID", "Ensembl_gene_identifier")
OUTPUT_COLUMN = "ensembl_gene_id"
LOG = logging.getLogger("add_ensembl_gene_ids")


def read_input_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read the input TSV and check that it has usable identifier columns."""
    with path.open(encoding="utf-8-sig", newline="") as input_file:
        reader = csv.DictReader(input_file, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Input TSV is empty or has no header")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("Input TSV contains duplicate column names")
        missing = REQUIRED_INPUT_COLUMNS - set(reader.fieldnames)
        if missing:
            raise ValueError(
                "Input TSV is missing required column(s): " + ", ".join(sorted(missing))
            )

        rows = []
        for line_number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Malformed TSV row at line {line_number}: column count differs from header")
            rows.append(row)
        return list(reader.fieldnames), rows


def gene_keys_from_rows(rows: list[dict[str, str]]) -> set[tuple[str, str]]:
    """Collect unique, usable (taxon ID, NCBI GeneID) keys from input rows."""
    return {
        (taxon_id, gene_id)
        for row in rows
        if (taxon_id := (row.get("taxon_id") or "").strip())
        and (gene_id := (row.get("ortholog_gene_id") or "").strip())
        and gene_id != "-"
    }


def load_gene_mappings(
    mapping_path: Path, wanted_keys: set[tuple[str, str]]
) -> dict[tuple[str, str], set[str]]:
    """Read only needed columns and keys from the compressed NCBI mapping TSV."""
    if not mapping_path.is_file():
        raise FileNotFoundError(
            f"NCBI mapping file not found: {mapping_path}\n"
            "Download it using the command in this experiment's README.md Downloads section."
        )

    if not wanted_keys:
        return {}

    mappings: dict[tuple[str, str], set[str]] = {}
    try:
        with gzip.open(mapping_path, "rt", encoding="utf-8-sig", newline="") as mapping_file:
            reader = csv.reader(mapping_file, delimiter="\t")
            try:
                header = next(reader)
            except StopIteration as error:
                raise ValueError("NCBI gene2ensembl.gz is empty") from error
            normalized_header = [name.lstrip("#").strip() for name in header]
            if len(normalized_header) != len(set(normalized_header)):
                raise ValueError("NCBI mapping file contains duplicate header names")
            missing = set(MAPPING_COLUMNS) - set(normalized_header)
            if missing:
                raise ValueError(
                    "NCBI mapping file is missing required column(s): "
                    + ", ".join(sorted(missing))
                )

            tax_index = normalized_header.index("tax_id")
            gene_index = normalized_header.index("GeneID")
            ensembl_index = normalized_header.index("Ensembl_gene_identifier")
            max_index = max(tax_index, gene_index, ensembl_index)
            wanted_taxa = {taxon_id for taxon_id, _ in wanted_keys}
            wanted_gene_ids = {gene_id for _, gene_id in wanted_keys}

            for line_number, fields in enumerate(reader, start=2):
                if len(fields) <= max_index:
                    raise ValueError(
                        f"Malformed NCBI mapping row at line {line_number}: "
                        "not enough tab-separated columns"
                    )
                taxon_id = fields[tax_index].strip()
                gene_id = fields[gene_index].strip()
                if taxon_id not in wanted_taxa or gene_id not in wanted_gene_ids:
                    continue
                key = (taxon_id, gene_id)
                if key not in wanted_keys:
                    continue
                ensembl_id = fields[ensembl_index].strip()
                if ensembl_id:
                    mappings.setdefault(key, set()).add(ensembl_id)
    except (OSError, EOFError) as error:
        raise ValueError(f"Could not read NCBI mapping file {mapping_path}: {error}") from error

    return mappings


def annotate_rows(
    rows: list[dict[str, str]], mappings: dict[tuple[str, str], set[str]]
) -> tuple[int, int, int]:
    """Set Ensembl IDs on each row and return single, multiple, and missing counts."""
    keys = gene_keys_from_rows(rows)
    single_count = multiple_count = missing_count = 0
    for key in keys:
        ensembl_ids = sorted(mappings.get(key, set()))
        if len(ensembl_ids) == 1:
            single_count += 1
        elif len(ensembl_ids) > 1:
            multiple_count += 1
            print(f"Multiple Ensembl mappings for taxon {key[0]}, GeneID {key[1]}: {';'.join(ensembl_ids)}")
        else:
            missing_count += 1
            print(f"No Ensembl mapping for taxon {key[0]}, GeneID {key[1]}")

    for row in rows:
        taxon_id = (row.get("taxon_id") or "").strip()
        gene_id = (row.get("ortholog_gene_id") or "").strip()
        ids = sorted(mappings.get((taxon_id, gene_id), set()))
        row[OUTPUT_COLUMN] = ";".join(ids)

    return single_count, multiple_count, missing_count


def automatic_output_path(input_path: Path) -> Path:
    """Return a sibling path with the input's ``01_`` prefix changed to ``02_``."""
    if not input_path.name.startswith("01_"):
        raise ValueError(
            f"Input filename must start with '01_' to generate the Step 2 output name: {input_path.name}"
        )
    return input_path.with_name("02_" + input_path.name[3:])


def write_output_table(
    output_path: Path, input_columns: list[str], rows: list[dict[str, str]]
) -> None:
    """Write all rows and input columns, appending the Ensembl ID column once."""
    columns = list(input_columns)
    if OUTPUT_COLUMN not in columns:
        columns.append(OUTPUT_COLUMN)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as output_file:
        writer = csv.DictWriter(
            output_file, fieldnames=columns, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def process_table(input_path: Path) -> Path:
    """Map all input orthologs locally and write the automatically named TSV."""
    input_columns, rows = read_input_table(input_path)
    output_path = automatic_output_path(input_path)
    keys = gene_keys_from_rows(rows)
    mappings = load_gene_mappings(MAPPING_FILE, keys)
    single_count, multiple_count, missing_count = annotate_rows(rows, mappings)
    write_output_table(output_path, input_columns, rows)

    print("\nSummary")
    print(f"Input rows: {len(rows)}")
    print(f"Unique NCBI genes processed: {len(keys)}")
    print(f"Genes with exactly one Ensembl mapping: {single_count}")
    print(f"Genes with multiple Ensembl mappings: {multiple_count}")
    print(f"Genes without an Ensembl mapping: {missing_count}")
    print(f"Output: {output_path}")
    return output_path


def main() -> int:
    """Parse the input path and run local Ensembl ID mapping."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_tsv", type=Path, help="Step 1 ortholog TSV")
    args = parser.parse_args()
    input_path = args.input_tsv.expanduser()
    try:
        process_table(input_path)
    except (OSError, ValueError, csv.Error) as error:
        LOG.error("Could not map Ensembl gene IDs: %s", error)
        return 1
    return 0


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    raise SystemExit(main())
