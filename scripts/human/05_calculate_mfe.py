#!/usr/bin/env python3
"""Calculate default ViennaRNA MFE structures for Step 3 start-codon windows."""

from __future__ import annotations

import csv
import argparse
import logging
import math
import multiprocessing
import os
import tempfile
from pathlib import Path
from typing import TextIO

import pandas as pd
import RNA
from Bio import SeqIO


ROOT = Path(__file__).resolve().parents[2]
INPUT_TSV = ROOT / "data" / "human" / "03_start_codon_windows.tsv"
INPUT_FASTA = ROOT / "data" / "human" / "03_start_codon_windows.fasta"
INTERMEDIATE_TSV = ROOT / "data" / "human" / "05_rnafold_results.tsv"
ERRORS_TSV = ROOT / "data" / "human" / "05_rnafold_errors.tsv"
FINAL_TSV = ROOT / "results" / "human" / "05_mane_start_codon_structure.tsv"
EXPECTED_REFERENCE_TRANSCRIPT = "NM_007278.2"
EXPECTED_REFERENCE_SEQUENCE = "TGTCCAGCCCGGTTCCCGGGAGGATGAAGTTCGTGTACAAAGAAGAGCAT"
EXPECTED_REFERENCE_GC = 54.0
EXPECTED_REFERENCE_MFE = -14.80
MFE_TOLERANCE = 0.01
PROGRESS_INTERVAL = 500
LOG = logging.getLogger("calculate_mfe")


def load_fasta_sequences(path: Path) -> dict[str, str]:
    """Load transcript sequences from the sequence-only FASTA source."""
    if not path.is_file():
        raise FileNotFoundError(f"Sequence FASTA not found: {path}")
    sequences: dict[str, str] = {}
    for record in SeqIO.parse(path, "fasta"):
        if record.id in sequences:
            raise ValueError(f"Duplicate transcript identifier in FASTA: {record.id}")
        sequences[record.id] = str(record.seq).upper()
    return sequences


def load_input(path: Path, fasta_path: Path) -> tuple[pd.DataFrame, dict[str, str]]:
    """Read metadata with pandas and sequence data from a TSV column or FASTA."""
    if not path.is_file():
        raise FileNotFoundError(f"Input TSV not found: {path}")
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"gene_name", "transcript", "GC"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Input TSV is missing required columns: {sorted(missing)}")
    if frame["transcript"].duplicated().any():
        duplicates = frame.loc[frame["transcript"].duplicated(keep=False), "transcript"].tolist()
        raise ValueError(f"Input contains duplicate transcript identifiers: {duplicates[:10]}")

    fasta_sequences: dict[str, str] = {}
    if "sequence" not in frame.columns or (frame["sequence"].str.strip() == "").any():
        fasta_sequences = load_fasta_sequences(fasta_path)
    return frame, fasta_sequences


def sequence_for_row(row: pd.Series, fasta_sequences: dict[str, str]) -> str:
    """Use an optional TSV sequence, otherwise look it up in the FASTA by ID."""
    sequence = row.get("sequence", "")
    if isinstance(sequence, str) and sequence.strip():
        return sequence.strip().upper()
    transcript = row["transcript"]
    if transcript not in fasta_sequences:
        raise ValueError(f"No sequence found in the TSV or FASTA for {transcript}")
    return fasta_sequences[transcript]


def fold_sequence(sequence_dna: str) -> tuple[str, float]:
    """Validate a 50-nt DNA window, convert T to U, and fold with defaults."""
    sequence_dna = sequence_dna.upper()
    if len(sequence_dna) != 50:
        raise ValueError(f"Expected 50 nt, received {len(sequence_dna)} nt")
    invalid = sorted(set(sequence_dna) - set("ACGT"))
    if invalid:
        raise ValueError(f"Invalid DNA symbols: {invalid}")
    sequence_rna = sequence_dna.replace("T", "U")
    structure, mfe = RNA.fold(sequence_rna)
    if len(structure) != len(sequence_rna):
        raise ValueError("ViennaRNA returned a structure with an unexpected length")
    return structure, float(mfe)


def process_transcript(task: tuple[str, str, str, str | None, str | None]) -> tuple[bool, dict[str, str]]:
    """Fold one transcript; kept at module scope for multiprocessing workers."""
    gene_name, transcript, gc_value, sequence_dna, input_error = task
    if input_error:
        error_type, description = input_error.split("\t", maxsplit=1)
        return False, {
            "gene_name": gene_name,
            "transcript": transcript,
            "error_type": error_type,
            "description": description,
        }
    try:
        structure, mfe = fold_sequence(sequence_dna or "")
        gc = float(gc_value)
        if not math.isfinite(gc) or not 0.0 <= gc <= 100.0:
            raise ValueError(f"Invalid GC percentage: {gc_value!r}")
        return True, {
            "gene_name": gene_name,
            "transcript": transcript,
            "structure": structure,
            "MFE": f"{mfe:.2f}",
            "_GC": f"{gc:.1f}",
            "_raw_MFE": str(mfe),
        }
    except Exception as error:
        return False, {
            "gene_name": gene_name,
            "transcript": transcript,
            "error_type": type(error).__name__,
            "description": str(error),
        }


def parse_args() -> argparse.Namespace:
    """Parse CLI options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--num_cpu",
        "--num-cpu",
        dest="num_cpu",
        type=int,
        default=1,
        help="Number of worker processes (default: 1; values above 1 enable multiprocessing)",
    )
    args = parser.parse_args()
    if args.num_cpu < 1:
        parser.error("--num_cpu must be at least 1")
    return args


def atomic_write_tsv(path: Path, rows: list[dict[str, str]], columns: tuple[str, ...]) -> None:
    """Write a TSV atomically with a fixed column order."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as handle:
        temporary_path = Path(handle.name)
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    temporary_path.replace(path)


def validate_gabarap(
    successful_rows: list[dict[str, str]], sequences_by_transcript: dict[str, str]
) -> None:
    """Compare the supplied GABARAP example with the computed result if present."""
    if EXPECTED_REFERENCE_TRANSCRIPT not in sequences_by_transcript:
        LOG.warning("GABARAP reference %s is absent from this input", EXPECTED_REFERENCE_TRANSCRIPT)
        return
    sequence = sequences_by_transcript[EXPECTED_REFERENCE_TRANSCRIPT]
    if sequence != EXPECTED_REFERENCE_SEQUENCE:
        LOG.warning("GABARAP reference sequence differs from the documented 50-nt example")
        return
    if (len(sequence) - sequence.count("A") - sequence.count("T")) / len(sequence) * 100 != EXPECTED_REFERENCE_GC:
        LOG.warning("GABARAP reference sequence does not have the expected 54%% GC")
    result = next(
        (row for row in successful_rows if row["transcript"] == EXPECTED_REFERENCE_TRANSCRIPT),
        None,
    )
    if result is None:
        LOG.warning("GABARAP reference sequence was not successfully folded")
        return
    mfe = float(result["_raw_MFE"])
    if abs(mfe - EXPECTED_REFERENCE_MFE) <= MFE_TOLERANCE:
        LOG.info("GABARAP MFE matches %.2f kcal/mol within %.2f", EXPECTED_REFERENCE_MFE, MFE_TOLERANCE)
    else:
        LOG.warning(
            "GABARAP MFE is %.2f kcal/mol; expected %.2f (tolerance %.2f). "
            "Check ViennaRNA version and folding settings; parameters were not changed.",
            mfe,
            EXPECTED_REFERENCE_MFE,
            MFE_TOLERANCE,
        )


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
    args = parse_args()
    version = getattr(RNA, "__version__", None) or getattr(RNA, "VERSION", "unknown")
    LOG.info("ViennaRNA Python API version: %s", version)
    LOG.info("Folding parameters: RNA.fold defaults; no constraints or custom model settings")
    LOG.info("Worker processes: %d", args.num_cpu)

    try:
        frame, fasta_sequences = load_input(INPUT_TSV, INPUT_FASTA)
        successful: list[dict[str, str]] = []
        errors: list[dict[str, str]] = []
        sequences_by_transcript: dict[str, str] = {}
        total = len(frame)

        tasks = []
        for _, row in frame.iterrows():
            transcript = row["transcript"].strip()
            gene_name = row["gene_name"].strip()
            try:
                sequence_dna = sequence_for_row(row, fasta_sequences)
                sequences_by_transcript[transcript] = sequence_dna
                tasks.append((gene_name, transcript, row["GC"], sequence_dna, None))
            except Exception as error:
                tasks.append(
                    (
                        gene_name,
                        transcript,
                        row["GC"],
                        None,
                        f"{type(error).__name__}\t{error}",
                    )
                )

        if args.num_cpu == 1:
            results = map(process_transcript, tasks)
            pool_context = None
        else:
            pool_context = multiprocessing.Pool(processes=args.num_cpu)
            chunk_size = max(1, min(50, len(tasks) // (args.num_cpu * 4) or 1))
            results = pool_context.imap(process_transcript, tasks, chunksize=chunk_size)

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
                if index % PROGRESS_INTERVAL == 0 or index == total:
                    LOG.info(
                        "Processed %d/%d transcripts; %d successes, %d errors",
                        index,
                        total,
                        len(successful),
                        len(errors),
                    )
        finally:
            if pool_context is not None:
                pool_context.close()
                pool_context.join()

        # Keep the public result tables limited to their specified metadata columns.
        intermediate_rows = [
            {key: row[key] for key in ("gene_name", "transcript", "structure", "MFE")}
            for row in successful
        ]
        final_rows = [
            {
                "gene_name": row["gene_name"],
                "transcript": row["transcript"],
                "GC": row["_GC"],
                "free_energy": row["MFE"],
            }
            for row in successful
        ]
        atomic_write_tsv(
            INTERMEDIATE_TSV,
            intermediate_rows,
            ("gene_name", "transcript", "structure", "MFE"),
        )
        atomic_write_tsv(
            FINAL_TSV,
            final_rows,
            ("gene_name", "transcript", "GC", "free_energy"),
        )
        atomic_write_tsv(
            ERRORS_TSV,
            errors,
            ("gene_name", "transcript", "error_type", "description"),
        )
        validate_gabarap(successful, sequences_by_transcript)
        LOG.info("Total input transcripts: %d", total)
        LOG.info("Successfully folded: %d", len(successful))
        LOG.info("Failed transcripts: %d", len(errors))
        LOG.info("Wrote intermediate results to %s", INTERMEDIATE_TSV)
        LOG.info("Wrote final results to %s", FINAL_TSV)
        LOG.info("Wrote error report to %s", ERRORS_TSV)
    except (OSError, ValueError, pd.errors.ParserError) as error:
        LOG.error("MFE calculation could not complete: %s", error)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
