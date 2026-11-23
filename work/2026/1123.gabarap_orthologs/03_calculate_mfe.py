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
INTERMEDIATE_COLUMNS = (
    "species_name",
    "species",
    "gene_symbol",
    "transcript",
    "structure",
    "MFE",
)
FINAL_COLUMNS = (
    "species_name",
    "species",
    "gene_symbol",
    "transcript",
    "GC",
    "free_energy",
)
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
) -> tuple[list[dict[str, str]], dict[str, str]]:
    """Read sequence metadata and load the FASTA only if TSV sequences are absent."""
    if not table_path.is_file():
        raise FileNotFoundError(f"Input TSV not found: {table_path}")
    with table_path.open(encoding="utf-8-sig", newline="") as input_file:
        reader = csv.DictReader(input_file, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Input TSV is empty or has no header")
        required = {"species_name", "species", "gene_symbol", "transcript", "GC"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(
                "Input TSV is missing required columns: "
                + ", ".join(sorted(missing))
            )
        rows = []
        identifiers = set()
        for line_number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(
                    f"Malformed TSV row at line {line_number}: "
                    "column count differs from header"
                )
            key = (row["species_name"], row["transcript"])
            if key in identifiers:
                raise ValueError(
                    f"Duplicate species/transcript pair in input TSV: {key}"
                )
            identifiers.add(key)
            rows.append(row)

    fasta_sequences: dict[str, str] = {}
    if any(not row.get("sequence", "").strip() for row in rows):
        fasta_sequences = load_fasta_sequences(fasta_path)
    return rows, fasta_sequences


def sequence_for_row(row: dict[str, str], fasta_sequences: dict[str, str]) -> str:
    """Use the TSV sequence or look it up in FASTA by species and transcript."""
    sequence = row.get("sequence", "").strip()
    if sequence:
        return sequence.upper()
    key = f"{row['species_name']}|{row['transcript']}"
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


def automatic_paths(gene_name: str) -> tuple[Path, Path, Path, Path, Path]:
    """Build default input and output paths from the required gene name."""
    safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", gene_name.strip()).strip("._-")
    if not safe_name:
        raise ValueError("Gene name must contain filename-safe characters")
    return (
        EXPERIMENT_DIR / f"02_{safe_name}_start_codon_windows.tsv",
        EXPERIMENT_DIR / f"02_{safe_name}_start_codon_windows.fasta",
        EXPERIMENT_DIR / f"03_{safe_name}_rnafold_results.tsv",
        EXPERIMENT_DIR / f"03_{safe_name}_start_codon_structure.tsv",
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
            and row["transcript"] == EXPECTED_HUMAN_TRANSCRIPT
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
            and row["transcript"] == EXPECTED_HUMAN_TRANSCRIPT
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
    parser.add_argument("--intermediate-output", type=Path, help="Override structure/MFE TSV")
    parser.add_argument("--final-output", type=Path, help="Override final energy TSV")
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
            default_intermediate,
            default_final,
            default_errors,
        ) = automatic_paths(args.gene_name)
        input_path = args.input or default_input
        fasta_path = args.fasta or default_fasta
        intermediate_path = args.intermediate_output or default_intermediate
        final_path = args.final_output or default_final
        errors_path = args.errors_output or default_errors

        import RNA

        version = getattr(RNA, "__version__", None) or getattr(RNA, "VERSION", "unknown")
        LOG.info("ViennaRNA Python API version: %s", version)
        LOG.info("Folding parameters: RNA.fold defaults; no constraints or custom model settings")
        LOG.info("Worker processes: %d", args.num_cpu)
        LOG.info("Input table: %s", input_path)

        rows, fasta_sequences = load_input(input_path, fasta_path)
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
                    row["transcript"],
                    row["GC"],
                    sequence,
                    input_error,
                )
            )

        successful: list[dict[str, str]] = []
        errors: list[dict[str, str]] = []
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

        intermediate_rows = [
            {key: row[key] for key in INTERMEDIATE_COLUMNS}
            for row in successful
        ]
        final_rows = [
            {
                "species_name": row["species_name"],
                "species": row["species"],
                "gene_symbol": row["gene_symbol"],
                "transcript": row["transcript"],
                "GC": row["_GC"],
                "free_energy": row["MFE"],
            }
            for row in successful
        ]
        atomic_write_tsv(intermediate_path, intermediate_rows, INTERMEDIATE_COLUMNS)
        atomic_write_tsv(final_path, final_rows, FINAL_COLUMNS)
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
        LOG.info("Wrote structure/MFE table: %s", intermediate_path)
        LOG.info("Wrote final energy table: %s", final_path)
    except (OSError, UnicodeError, csv.Error, ValueError, ImportError) as error:
        LOG.error("MFE calculation could not complete: %s", error)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
