#!/usr/bin/env python3
"""Add Ensembl canonical transcript data from local preprocessed TSV tables."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


EXPERIMENT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = EXPERIMENT_DIR.parents[2]
OUTPUT_COLUMNS = (
    "ensembl_gene",
    "ensembl_canonical_transcript",
    "ensembl_protein",
    "cds_start",
    "cds_end",
    "cds_length",
)
CANONICAL_COLUMNS = {
    "gene_id",
    "gene_name",
    "gene_biotype",
    "transcript_id",
    "transcript_biotype",
    "protein_id",
    "cds_start",
    "cds_end",
    "cds_length",
}
PROGRESS_INTERVAL = 100_000

sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "mammals"))
from constants import SPECIES_DICT  # noqa: E402


@dataclass(frozen=True)
class MatchResult:
    """Ensembl identifiers and CDS fields resolved for one gene symbol."""

    values: tuple[str, str, str, str, str, str] = ("-",) * 6
    status: str = "missing_gene"
    reason: str = ""


def read_input_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read and validate the ortholog input TSV."""
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


def canonical_table_for_species(species_name: str, table_dir: Path) -> Path:
    """Find the single local canonical table for a SPECIES_DICT key."""
    if species_name not in SPECIES_DICT:
        raise ValueError(
            f"Species key {species_name!r} is not present in SPECIES_DICT"
        )
    scientific_name = str(SPECIES_DICT[species_name].get("species", "")).strip()
    if not scientific_name:
        raise ValueError(f"Species key {species_name!r} has no scientific name")
    prefix = scientific_name.replace(" ", "_")
    prefix = prefix[0].upper() + prefix[1:]
    matches = sorted(table_dir.glob(f"{prefix}.*.canonical_transcripts.tsv"))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected one canonical transcript table for {species_name!r}, "
            f"matching {prefix}.*.canonical_transcripts.tsv in {table_dir}; "
            f"found {len(matches)}"
        )
    return matches[0]


def read_canonical_table(
    table_path: Path, requested_symbols: set[str]
) -> dict[str, dict[str, dict[str, dict[str, str]]]]:
    """Read matching genes and transcripts from one canonical transcript TSV."""
    matches: dict[str, dict[str, dict[str, dict[str, str]]]] = defaultdict(
        lambda: defaultdict(dict)
    )
    line_count = 0
    print(f"Reading canonical transcript table: {table_path.name}", flush=True)
    with table_path.open(encoding="utf-8-sig", newline="") as table_file:
        reader = csv.DictReader(table_file, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Canonical transcript table has no header: {table_path}")
        missing = CANONICAL_COLUMNS - set(reader.fieldnames)
        if missing:
            raise ValueError(
                f"Canonical transcript table {table_path} is missing columns: "
                + ", ".join(sorted(missing))
            )
        for line_count, row in enumerate(reader, start=1):
            if line_count % PROGRESS_INTERVAL == 0:
                print(
                    f"  Processed {line_count:,} transcript rows; "
                    f"{sum(len(genes) for genes in matches.values()):,} "
                    "matching gene(s).",
                    flush=True,
                )
            if None in row or any(value is None for value in row.values()):
                raise ValueError(
                    f"Malformed row {line_count + 1} in {table_path}: "
                    "column count differs from header"
                )
            symbol = row["gene_name"]
            if symbol not in requested_symbols:
                continue
            if (
                row["gene_biotype"] != "protein_coding"
                or row["transcript_biotype"] != "protein_coding"
            ):
                continue
            gene_id = row["gene_id"]
            transcript_id = row["transcript_id"]
            if not gene_id or not transcript_id:
                raise ValueError(
                    f"Missing gene_id or transcript_id in {table_path}, "
                    f"data row {line_count}"
                )
            matches[symbol][gene_id][transcript_id] = row

    print(
        f"  Read {line_count:,} transcript rows; "
        f"{sum(len(genes) for genes in matches.values()):,} matching gene(s), "
        f"{sum(len(transcripts) for genes in matches.values() for transcripts in genes.values()):,} "
        "canonical protein-coding transcript(s).",
        flush=True,
    )
    return matches


def resolve_symbol(
    symbol: str,
    matches: dict[str, dict[str, dict[str, dict[str, str]]]],
) -> MatchResult:
    """Resolve exact gene-name matches without choosing ambiguous IDs."""
    genes = matches.get(symbol, {})
    if not symbol or symbol == "-" or not genes:
        return MatchResult(reason="no matching gene in the canonical transcript table")
    if len(genes) > 1:
        return MatchResult(
            status="ambiguous_gene",
            reason=f"multiple gene IDs have gene_name {symbol!r}",
        )

    gene_id, transcripts = next(iter(genes.items()))
    if not transcripts:
        return MatchResult(
            status="missing_transcript",
            reason="no canonical protein-coding transcript in the table",
        )
    if len(transcripts) > 1:
        return MatchResult(
            values=(gene_id, "-", "-", "-", "-", "-"),
            status="ambiguous_transcript",
            reason="multiple canonical protein-coding transcripts match the gene",
        )

    transcript_id, row = next(iter(transcripts.items()))
    protein_id = row["protein_id"] or "-"
    cds_values = [row[column] for column in ("cds_start", "cds_end", "cds_length")]
    cds_status = "matched"
    cds_reason = ""
    if any(value in {"", "-"} for value in cds_values):
        cds_values = ["-", "-", "-"]
        cds_status = "missing_cds"
        cds_reason = "CDS coordinates are missing"
    else:
        try:
            cds_start, cds_end, cds_length = (int(value) for value in cds_values)
        except ValueError:
            cds_values = ["-", "-", "-"]
            cds_status = "invalid_cds"
            cds_reason = "CDS coordinates are not integers"
        else:
            if (
                cds_start < 1
                or cds_end < cds_start
                or cds_length != cds_end - cds_start + 1
            ):
                cds_values = ["-", "-", "-"]
                cds_status = "invalid_cds"
                cds_reason = "CDS coordinate range or length is invalid"

    if cds_status != "matched":
        status = cds_status
    elif protein_id == "-":
        status = "missing_protein"
    else:
        status = "complete"
    protein_reason = "missing protein_id" if protein_id == "-" else ""
    reasons = [reason for reason in (cds_reason, protein_reason) if reason]
    return MatchResult(
        values=(gene_id, transcript_id, protein_id, *cds_values),
        status=status,
        reason="; ".join(reasons),
    )


def write_output(
    path: Path, input_columns: list[str], rows: list[dict[str, str]]
) -> None:
    """Write all input records and append Ensembl annotations."""
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
        writer.writerows(
            {column: row.get(column) or "-" for column in columns}
            for row in rows
        )


def automatic_paths(gene_name: str) -> tuple[Path, Path]:
    """Build the default input and output filenames from a gene name."""
    safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", gene_name.strip()).strip("._-")
    if not safe_name:
        raise ValueError("Gene name must contain filename-safe characters")
    return (
        EXPERIMENT_DIR / f"00_{safe_name}_orthologs.tsv",
        EXPERIMENT_DIR / f"01_{safe_name}_orthologs.tsv",
    )


def process_table(
    input_path: Path, output_path: Path, canonical_dir: Path
) -> int:
    """Annotate each input row using its species' local canonical table."""
    print(f"Reading input table: {input_path}", flush=True)
    input_columns, rows = read_input_table(input_path)
    print(f"Input rows: {len(rows):,}", flush=True)

    rows_by_species: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        species_name = row["species_name"]
        if species_name not in SPECIES_DICT:
            raise ValueError(
                f"Species key {species_name!r} is not present in SPECIES_DICT"
            )
        rows_by_species[species_name].append(row)

    overall = Counter()
    for species_name, species_rows in rows_by_species.items():
        table_path = canonical_table_for_species(species_name, canonical_dir)
        symbols = {
            row["ortholog_gene_symbol"]
            for row in species_rows
            if row["ortholog_gene_symbol"] not in {"", "-"}
        }
        print(
            f"\nProcessing {species_name}: {len(species_rows)} row(s), "
            f"{len(symbols)} unique symbol(s).",
            flush=True,
        )
        matches = read_canonical_table(table_path, symbols)
        species_counts = Counter()
        results = {symbol: resolve_symbol(symbol, matches) for symbol in symbols}

        for row in species_rows:
            symbol = row["ortholog_gene_symbol"]
            result = results.get(symbol) or resolve_symbol(symbol, matches)
            row.update(zip(OUTPUT_COLUMNS, result.values))
            species_counts[result.status] += 1
            species_counts["rows"] += 1
            if result.reason:
                print(
                    f"  {symbol!r}: {result.status} ({result.reason}).",
                    flush=True,
                )

        print(f"Species summary ({species_name}):", flush=True)
        print(f"  Rows: {species_counts['rows']}", flush=True)
        print(f"  Complete matches: {species_counts['complete']}", flush=True)
        print(f"  Missing genes: {species_counts['missing_gene']}", flush=True)
        print(f"  Ambiguous gene matches: {species_counts['ambiguous_gene']}", flush=True)
        print(
            f"  Ambiguous canonical transcripts: "
            f"{species_counts['ambiguous_transcript']}",
            flush=True,
        )
        print(f"  Missing CDS: {species_counts['missing_cds']}", flush=True)
        print(f"  Invalid CDS: {species_counts['invalid_cds']}", flush=True)
        print(f"  Missing protein IDs: {species_counts['missing_protein']}", flush=True)
        overall.update(species_counts)

    print(f"\nWriting output table: {output_path}", flush=True)
    write_output(output_path, input_columns, rows)
    print("Final summary across species:", flush=True)
    print(f"  Input rows: {overall['rows']}", flush=True)
    print(f"  Complete matches: {overall['complete']}", flush=True)
    print(f"  Missing genes: {overall['missing_gene']}", flush=True)
    print(f"  Ambiguous gene matches: {overall['ambiguous_gene']}", flush=True)
    print(
        f"  Ambiguous canonical transcripts: {overall['ambiguous_transcript']}",
        flush=True,
    )
    print(f"  Missing CDS: {overall['missing_cds']}", flush=True)
    print(f"  Invalid CDS: {overall['invalid_cds']}", flush=True)
    print(f"  Missing protein IDs: {overall['missing_protein']}", flush=True)
    print(f"  Output: {output_path}", flush=True)
    return 0


def main() -> int:
    """Parse arguments and add local Ensembl canonical transcript annotations."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "gene_name",
        help="Gene name used to select the default input and output filenames",
    )
    parser.add_argument("--input", type=Path, help="Override the input TSV path")
    parser.add_argument("--output", type=Path, help="Override the output TSV path")
    parser.add_argument(
        "--canonical-dir",
        type=Path,
        default=EXPERIMENT_DIR,
        help="Folder containing species canonical_transcripts.tsv tables",
    )
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        return 2
    args = parser.parse_args()
    try:
        default_input, default_output = automatic_paths(args.gene_name)
        return process_table(
            args.input or default_input,
            args.output or default_output,
            args.canonical_dir,
        )
    except (OSError, UnicodeError, csv.Error, ValueError) as error:
        parser.error(f"Could not add canonical transcript annotations: {error}")


if __name__ == "__main__":
    raise SystemExit(main())
