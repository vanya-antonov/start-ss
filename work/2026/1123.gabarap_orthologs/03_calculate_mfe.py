#!/usr/bin/env python3
"""Calculate default ViennaRNA minimum-free-energy folds for start-codon windows."""

from __future__ import annotations

import argparse
import csv
import logging
import math
import multiprocessing
import os
import re
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import TextIO


EXPERIMENT_DIR = Path(__file__).resolve().parent
EXPECTED_HUMAN_TRANSCRIPT = "ENST00000302386.10"
EXPECTED_HUMAN_SEQUENCE = "TGTCCAGCCCGGTTCCCGGGAGGATGAAGTTCGTGTACAAAGAAGAGCAT"
EXPECTED_HUMAN_GC = 54.0
EXPECTED_HUMAN_MFE = -14.80
MFE_TOLERANCE = 0.01
PROGRESS_INTERVAL = 500
LOG = logging.getLogger("calculate_gabarap_mfe")
ERROR_COLUMNS = (
    "species_name",
    "species",
    "gene_symbol",
    "transcript",
    "error_type",
    "description",
)


def fold_sequence(sequence_dna: str) -> tuple[str, float]:
    """Validate a 50-nt DNA window and fold its RNA equivalent with defaults."""
    sequence_dna = sequence_dna.upper()
    if len(sequence_dna) != 50:
        raise ValueError(f"Expected 50 nt, received {len(sequence_dna)} nt")
    invalid = sorted(set(sequence_dna) - set("ACGT"))
    if invalid:
        raise ValueError(f"Invalid DNA symbols: {invalid}")
    sequence_rna = sequence_dna.replace("T", "U")
    import RNA

    structure, mfe = RNA.fold(sequence_rna)
    if len(structure) != len(sequence_rna):
        raise ValueError("ViennaRNA returned a structure with an unexpected length")
    return structure, float(mfe)


def load_fasta_sequences(path: Path) -> dict[str, str]:
    """Read sequences keyed by the FASTA record ID."""
    if not path.is_file():
        raise FileNotFoundError(f"Sequence FASTA not found: {path}")
    sequences: dict[str, str] = {}
    current_id = ""
    parts: list[str] = []

    def save_record() -> None:
        if not current_id:
            return
        if current_id in sequences:
            raise ValueError(f"Duplicate sequence identifier in FASTA: {current_id}")
        sequences[current_id] = "".join(parts).upper()

    with path.open(encoding="utf-8") as fasta_file:
        for line_number, line in enumerate(fasta_file, start=1):
            if line.startswith(">"):
                save_record()
                header = line[1:].strip()
                current_id = header.split(maxsplit=1)[0] if header else ""
                parts = []
                if not current_id:
                    raise ValueError(f"Empty FASTA header at line {line_number}")
            elif current_id:
                parts.append("".join(line.split()))
            elif line.strip():
                raise ValueError(
                    f"Sequence data before the first FASTA header at line {line_number}"
                )
    save_record()
    return sequences


def load_input(
    table_path: Path, fasta_path: Path
) -> tuple[list[str], list[dict[str, str]], dict[str, str]]:
    """Read sequence metadata and load the FASTA only if TSV sequences are absent."""
    if not table_path.is_file():
        raise FileNotFoundError(f"Input TSV not found: {table_path}")
    with table_path.open(encoding="utf-8-sig", newline="") as input_file:
        reader = csv.DictReader(input_file, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Input TSV is empty or has no header")
        required = {"species_name", "species", "gene_symbol", "GC"}
        missing = required - set(reader.fieldnames)
        if not ({"transcript", "ensembl_canonical_transcript"} & set(reader.fieldnames)):
            missing.add("transcript or ensembl_canonical_transcript")
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

    fasta_sequences: dict[str, str] = {}
    if any(row.get("sequence", "").strip() in {"", "-"} for row in rows):
        fasta_sequences = load_fasta_sequences(fasta_path)
    return list(reader.fieldnames), rows, fasta_sequences


def transcript_for_row(row: dict[str, str]) -> str:
    """Return the transcript ID column used by the input table."""
    return row.get("transcript") or row.get("ensembl_canonical_transcript", "")


def sequence_for_row(row: dict[str, str], fasta_sequences: dict[str, str]) -> str:
    """Use the TSV sequence or look it up in FASTA by species and transcript."""
    sequence = row.get("sequence", "").strip()
    if sequence and sequence != "-":
        return sequence.upper()
    transcript = transcript_for_row(row)
    key = f"{row['species_name']}|{transcript}"
    if key not in fasta_sequences:
        raise ValueError(f"No sequence found in TSV or FASTA for {key}")
    return fasta_sequences[key]


def process_transcript(
    task: tuple[str, str, str, str, str, str | None, str | None]
) -> tuple[bool, dict[str, str]]:
    """Fold one sequence; this module-level function can run in worker processes."""
    species_name, species, gene_symbol, transcript, gc_value, sequence_dna, input_error = task
    base = {
        "species_name": species_name,
        "species": species,
        "gene_symbol": gene_symbol,
        "transcript": transcript,
    }
    if input_error:
        error_type, description = input_error.split("\t", maxsplit=1)
        return False, {
            **base,
            "error_type": error_type,
            "description": description,
        }

    try:
        structure, mfe = fold_sequence(sequence_dna or "")
        gc = float(gc_value)
        if not math.isfinite(gc) or not 0.0 <= gc <= 100.0:
            raise ValueError(f"Invalid GC percentage: {gc_value!r}")
        return True, {
            **base,
            "structure": structure,
            "MFE": f"{mfe:.2f}",
            "_GC": f"{gc:.1f}",
            "_raw_MFE": str(mfe),
        }
    except Exception as error:
        return False, {
            **base,
            "error_type": type(error).__name__,
            "description": str(error),
        }


def atomic_write_tsv(
    path: Path, rows: list[dict[str, str]], columns: tuple[str, ...]
) -> None:
    """Write a TSV atomically with a fixed column order."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as handle:
        temporary_path = Path(handle.name)
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    temporary_path.replace(path)


def automatic_paths(gene_name: str) -> tuple[Path, Path, Path, Path]:
    """Build default input and output paths from the required gene name."""
    safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", gene_name.strip()).strip("._-")
    if not safe_name:
        raise ValueError("Gene name must contain filename-safe characters")
    return (
        EXPERIMENT_DIR / f"02_{safe_name}_orthologs.tsv",
        EXPERIMENT_DIR / f"02_{safe_name}_start_codon_windows.fasta",
        EXPERIMENT_DIR / f"03_{safe_name}_orthologs.tsv",
        EXPERIMENT_DIR / f"03_{safe_name}_rnafold_errors.tsv",
    )


def validate_human_reference(
    successful_rows: list[dict[str, str]], rows: list[dict[str, str]]
) -> None:
    """Check the known human GABARAP sequence and its default MFE when present."""
    reference = next(
        (
            row
            for row in rows
            if row["species_name"] == "human"
            and transcript_for_row(row) == EXPECTED_HUMAN_TRANSCRIPT
        ),
        None,
    )
    if reference is None:
        LOG.warning("Human reference transcript %s is absent from input", EXPECTED_HUMAN_TRANSCRIPT)
        return
    sequence = reference.get("sequence", "").strip().upper()
    if sequence and sequence != EXPECTED_HUMAN_SEQUENCE:
        LOG.warning("Human reference sequence differs from the documented 50-nt example")
        return
    result = next(
        (
            row
            for row in successful_rows
            if row["species_name"] == "human"
            and transcript_for_row(row) == EXPECTED_HUMAN_TRANSCRIPT
        ),
        None,
    )
    if result is None:
        LOG.warning("Human reference sequence was not successfully folded")
        return
    mfe = float(result["_raw_MFE"])
    if abs(mfe - EXPECTED_HUMAN_MFE) <= MFE_TOLERANCE:
        LOG.info(
            "Human GABARAP MFE matches %.2f kcal/mol within %.2f",
            EXPECTED_HUMAN_MFE,
            MFE_TOLERANCE,
        )
    else:
        LOG.warning(
            "Human GABARAP MFE is %.2f kcal/mol; expected %.2f (tolerance %.2f). "
            "Check ViennaRNA version and folding settings; parameters were not changed.",
            mfe,
            EXPECTED_HUMAN_MFE,
            MFE_TOLERANCE,
        )


def report_gc_distribution(rows: list[dict[str, str]]) -> None:
    """Log compact summary statistics for the GC percentages."""
    values = [float(row["_GC"]) for row in rows]
    if not values:
        LOG.info("GC content distribution: no successful sequences")
        return
    values.sort()
    midpoint = len(values) // 2
    median = (
        values[midpoint]
        if len(values) % 2
        else (values[midpoint - 1] + values[midpoint]) / 2
    )
    LOG.info(
        "GC%% distribution: min=%.1f, median=%.1f, mean=%.1f, max=%.1f",
        values[0],
        median,
        sum(values) / len(values),
        values[-1],
    )
    for lower, upper in ((0, 40), (40, 50), (50, 60), (60, 101)):
        count = sum(lower <= value < upper for value in values)
        label = f"{lower}-<{upper}%" if upper != 101 else f"{lower}%+"
        LOG.info("GC%% %s: %d", label, count)


def parse_args() -> argparse.Namespace:
    """Parse the required gene name and optional paths/worker count."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "gene_name",
        help="Gene name used to select default input and output filenames",
    )
    parser.add_argument("--input", type=Path, help="Override the input results TSV")
    parser.add_argument("--fasta", type=Path, help="Override the sequence FASTA path")
    parser.add_argument(
        "--num-cpu",
        "--num_cpu",
        dest="num_cpu",
        type=int,
        default=1,
        help="Number of worker processes (default: 1)",
    )
    parser.add_argument("--output", type=Path, help="Override output TSV")
    parser.add_argument("--errors-output", type=Path, help="Override error report TSV")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        raise SystemExit(2)
    args = parser.parse_args()
    if args.num_cpu < 1:
        parser.error("--num-cpu must be at least 1")
    return args


def main() -> int:
    """Calculate default ViennaRNA MFE structures and write result tables."""
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s"
    )
    args = parse_args()
    try:
        (
            default_input,
            default_fasta,
            default_output,
            default_errors,
        ) = automatic_paths(args.gene_name)
        input_path = args.input or default_input
        fasta_path = args.fasta or default_fasta
        output_path = args.output or default_output
        errors_path = args.errors_output or default_errors

        import RNA

        version = getattr(RNA, "__version__", None) or getattr(RNA, "VERSION", "unknown")
        LOG.info("ViennaRNA Python API version: %s", version)
        LOG.info("Folding parameters: RNA.fold defaults; no constraints or custom model settings")
        LOG.info("Worker processes: %d", args.num_cpu)
        LOG.info("Input table: %s", input_path)

        input_columns, rows, fasta_sequences = load_input(input_path, fasta_path)
        tasks = []
        for row in rows:
            try:
                sequence = sequence_for_row(row, fasta_sequences)
                input_error = None
            except Exception as error:
                sequence = None
                input_error = f"{type(error).__name__}\t{error}"
            tasks.append(
                (
                    row["species_name"],
                    row["species"],
                    row["gene_symbol"],
                    transcript_for_row(row),
                    row["GC"],
                    sequence,
                    input_error,
                )
            )

        successful: list[dict[str, str]] = []
        errors: list[dict[str, str]] = []
        output_rows = [
            {**row, "structure": "-", "MFE": "-"} for row in rows
        ]
        pool = None
        if args.num_cpu == 1:
            results = map(process_transcript, tasks)
        else:
            pool = multiprocessing.Pool(processes=args.num_cpu)
            chunk_size = max(1, min(50, len(tasks) // (args.num_cpu * 4) or 1))
            results = pool.imap(process_transcript, tasks, chunksize=chunk_size)

        try:
            for index, (is_success, result) in enumerate(results, start=1):
                if is_success:
                    successful.append(result)
                    output_rows[index - 1]["structure"] = result["structure"]
                    output_rows[index - 1]["MFE"] = result["MFE"]
                else:
                    errors.append(result)
                    LOG.error(
                        "Failed %s (%s): %s",
                        result["transcript"],
                        result["error_type"],
                        result["description"],
                    )
                if index % PROGRESS_INTERVAL == 0 or index == len(tasks):
                    LOG.info(
                        "Processed %d/%d transcripts; %d successes, %d errors",
                        index,
                        len(tasks),
                        len(successful),
                        len(errors),
                    )
        finally:
            if pool is not None:
                pool.close()
                pool.join()

        output_columns = list(input_columns)
        for column in ("structure", "MFE"):
            if column not in output_columns:
                output_columns.append(column)
        atomic_write_tsv(output_path, output_rows, tuple(output_columns))
        if errors:
            atomic_write_tsv(errors_path, errors, ERROR_COLUMNS)
            LOG.info("Wrote error report: %s", errors_path)
        else:
            errors_path.unlink(missing_ok=True)
            LOG.info("No folding errors; error report not created")

        validate_human_reference(successful, rows)
        report_gc_distribution(successful)
        LOG.info("Total input windows: %d", len(tasks))
        LOG.info("Successfully folded: %d", len(successful))
        LOG.info("Failed windows: %d", len(errors))
        LOG.info("Wrote structure/MFE table: %s", output_path)
    except (OSError, UnicodeError, csv.Error, ValueError, ImportError) as error:
        LOG.error("MFE calculation could not complete: %s", error)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
