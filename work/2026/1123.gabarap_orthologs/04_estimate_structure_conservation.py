#!/usr/bin/env python3
"""Estimate dot-bracket structure conservation against the human transcript."""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import tempfile
from pathlib import Path
EXPERIMENT_DIR = Path(__file__).resolve().parent


def read_structures(path: Path) -> list[dict[str, str]]:
    """Read rows and require the fields needed for comparison."""
    if not path.is_file():
        raise FileNotFoundError(f"Input TSV not found: {path}")
    with path.open(encoding="utf-8-sig", newline="") as input_file:
        reader = csv.DictReader(input_file, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Input TSV is empty or has no header")
        required = {"species_name", "species", "structure"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(
                "Input TSV is missing required columns: "
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
    return rows


def is_valid_structure(structure: str) -> bool:
    """Return whether a value is a nonempty, valid dot-bracket structure."""
    structure = structure.strip()
    if not structure or structure == "-" or set(structure) - set(".()"):
        return False
    depth = 0
    for symbol in structure:
        if symbol == "(":
            depth += 1
        elif symbol == ")":
            depth -= 1
            if depth < 0:
                return False
    return depth == 0


def transcript_for_row(row: dict[str, str]) -> str:
    """Return the transcript identifier from either supported column name."""
    return row.get("transcript") or row.get("ensembl_canonical_transcript", "-")


def atomic_write_tsv(
    path: Path, rows: list[dict[str, str]], columns: tuple[str, ...] | list[str]
) -> None:
    """Write a TSV atomically."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as output_file:
        temporary_path = Path(output_file.name)
        writer = csv.DictWriter(
            output_file, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
        output_file.flush()
        os.fsync(output_file.fileno())
    temporary_path.replace(path)


def automatic_paths(gene_name: str) -> tuple[Path, Path, Path]:
    """Build default input and output paths from the required gene name."""
    safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", gene_name.strip()).strip("._-")
    if not safe_name:
        raise ValueError("Gene name must contain filename-safe characters")
    return (
        EXPERIMENT_DIR / f"03_{safe_name}_orthologs.tsv",
        EXPERIMENT_DIR / f"04_{safe_name}_structure_conservation.tsv",
        EXPERIMENT_DIR / f"04_{safe_name}_orthologs.tsv",
    )


def compare_to_human(
    rows: list[dict[str, str]], max_mismatches: int
) -> dict[str, str]:
    """Compare equal-length nonhuman structures to the single human structure."""
    human_rows = [row for row in rows if row["species_name"].strip().lower() == "human"]
    if len(human_rows) != 1:
        raise ValueError(f"Expected exactly one human row; found {len(human_rows)}")
    human = human_rows[0]
    human_structure = human["structure"].strip()
    if not is_valid_structure(human_structure):
        raise ValueError("The human row does not contain a valid dot-bracket structure")

    ortholog_rows = [row for row in rows if row is not human]
    valid_orthologs = [
        row for row in ortholog_rows if is_valid_structure(row["structure"])
    ]
    if not valid_orthologs:
        raise ValueError("No nonhuman rows have valid dot-bracket structures")

    expected_length = len(human_structure)
    length_mismatches = [
        (row["species_name"], len(row["structure"].strip()))
        for row in valid_orthologs
        if len(row["structure"].strip()) != expected_length
    ]
    if length_mismatches:
        details = ", ".join(f"{species}: {length}" for species, length in length_mismatches)
        raise ValueError(
            f"All valid structures must have length {expected_length}; "
            f"mismatch(es): {details}"
        )

    conserved = 0
    for row in rows:
        row["num_structure_mismatches_with_human"] = "-"
        row["structure_match_human"] = "-"

    for row in valid_orthologs:
        structure = row["structure"].strip()
        mismatches = sum(a != b for a, b in zip(human_structure, structure))
        row["num_structure_mismatches_with_human"] = str(mismatches)
        row["structure_match_human"] = "yes" if mismatches <= max_mismatches else "no"
        if mismatches <= max_mismatches:
            conserved += 1

    percentage = 100.0 * conserved / len(valid_orthologs)
    human_transcript = transcript_for_row(human)
    return {
        "reference_species": human["species"],
        "reference_transcript": human_transcript or "-",
        "reference_structure_length": str(expected_length),
        "num_smiles_mismatches": str(max_mismatches),
        "num_orthologs": str(len(valid_orthologs)),
        "num_conserved_orthologs": str(conserved),
        "percent_conserved": f"{percentage:.1f}",
        "nonhuman_rows_excluded_missing_or_invalid_structure": str(
            len(ortholog_rows) - len(valid_orthologs)
        ),
    }


def parse_args() -> argparse.Namespace:
    """Parse the required input file and mismatch threshold."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "gene_name",
        help="Gene name used to select default input and output filenames",
    )
    parser.add_argument(
        "--num_smiles_mismatches",
        type=int,
        default=2,
        help="Maximum differing dot-bracket characters for a conserved structure",
    )
    parser.add_argument("--input", type=Path, help="Override the Step 3 input TSV")
    parser.add_argument("--output", type=Path, help="Override the summary TSV path")
    parser.add_argument(
        "--orthologs-output",
        type=Path,
        help="Override the annotated ortholog TSV path",
    )
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        raise SystemExit(2)
    args = parser.parse_args()
    if args.num_smiles_mismatches < 0:
        parser.error("--num_smiles_mismatches must be zero or greater")
    return args


def main() -> int:
    """Estimate and report structure conservation relative to human."""
    args = parse_args()
    try:
        default_input, default_output, default_orthologs_output = automatic_paths(
            args.gene_name
        )
        input_path = args.input or default_input
        output_path = args.output or default_output
        orthologs_output_path = args.orthologs_output or default_orthologs_output
        print(f"Reading Step 3 structures: {input_path}", flush=True)
        rows = read_structures(input_path)
        print(f"Input rows: {len(rows):,}", flush=True)
        summary = compare_to_human(rows, args.num_smiles_mismatches)
        summary_rows = [
            {"statistic": name, "value": value} for name, value in summary.items()
        ]
        atomic_write_tsv(output_path, summary_rows, ("statistic", "value"))
        ortholog_columns = list(rows[0])
        for column in (
            "num_structure_mismatches_with_human",
            "structure_match_human",
        ):
            if column not in ortholog_columns:
                ortholog_columns.append(column)
        atomic_write_tsv(orthologs_output_path, rows, ortholog_columns)
        print(
            f"Compared {summary['num_orthologs']} valid nonhuman structure(s); "
            f"{summary['num_conserved_orthologs']} conserved "
            f"({summary['percent_conserved']}%).",
            flush=True,
        )
        print(
            "Excluded nonhuman rows with missing or invalid structures: "
            f"{summary['nonhuman_rows_excluded_missing_or_invalid_structure']}",
            flush=True,
        )
        print(f"Summary TSV: {output_path}", flush=True)
        print(f"Annotated orthologs TSV: {orthologs_output_path}", flush=True)
    except (OSError, UnicodeError, csv.Error, ValueError) as error:
        print(f"Structure conservation could not be calculated: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
