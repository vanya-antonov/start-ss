#!/usr/bin/env python3
"""Extract validated 50-nt windows centered on annotated CDS start codons."""

from __future__ import annotations

import csv
import gzip
import logging
import os
import tempfile
from collections import Counter
from pathlib import Path
from statistics import mean, median
from typing import Callable, TextIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord


DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "human"
INPUT_FILE = DATA_DIR / "02_mane_transcripts.gb.gz"
FASTA_FILE = DATA_DIR / "03_start_codon_windows.fasta"
RESULTS_FILE = DATA_DIR / "03_start_codon_windows.tsv"
ERRORS_FILE = DATA_DIR / "03_extraction_errors.tsv"
UPSTREAM_LENGTH = 23
START_CODON_LENGTH = 3
DOWNSTREAM_LENGTH = 24
WINDOW_LENGTH = UPSTREAM_LENGTH + START_CODON_LENGTH + DOWNSTREAM_LENGTH
VALID_BASES = frozenset("ACGT")
REFERENCE_TRANSCRIPT = "NM_007278.2"
REFERENCE_SEQUENCE = "TGTCCAGCCCGGTTCCCGGGAGGATGAAGTTCGTGTACAAAGAAGAGCAT"
LOG = logging.getLogger("extract_start_sequences")


class ExtractionError(Exception):
    """An expected per-transcript extraction failure with a reportable type."""

    def __init__(self, error_type: str, description: str) -> None:
        super().__init__(description)
        self.error_type = error_type


def gene_name_for_record(record: SeqRecord, cds_features: list[SeqFeature]) -> str:
    """Get a gene symbol from CDS/gene feature qualifiers when available."""
    names: list[str] = []
    for feature in cds_features:
        names.extend(feature.qualifiers.get("gene", []))
    if not names:
        for feature in record.features:
            if feature.type == "gene":
                names.extend(feature.qualifiers.get("gene", []))
    unique_names = list(dict.fromkeys(name.strip() for name in names if name.strip()))
    return ";".join(unique_names)


def transcript_accession(record: SeqRecord) -> str:
    """Return the versioned accession represented by the GenBank record ID."""
    return record.id.strip()


def extract_window(
    record: SeqRecord, cds_features: list[SeqFeature]
) -> tuple[str, str, str, float]:
    """Extract a validated start-codon window without changing transcript strand."""
    transcript = transcript_accession(record)
    gene_name = gene_name_for_record(record, cds_features)
    if not cds_features:
        raise ExtractionError("missing_cds", "No CDS feature is annotated")
    if len(cds_features) > 1:
        raise ExtractionError(
            "multiple_cds", f"Found {len(cds_features)} CDS features; cannot choose one unambiguously"
        )

    sequence = str(record.seq).upper()
    # RefSeq transcript records are already in transcript (5′ to 3′) orientation.
    # CDS location coordinates therefore index the stored sequence directly.
    cds_start = int(cds_features[0].location.start)
    if cds_start < UPSTREAM_LENGTH:
        raise ExtractionError(
            "insufficient_upstream_sequence",
            f"CDS starts at position {cds_start}; need {UPSTREAM_LENGTH} upstream nucleotides",
        )
    if cds_start + START_CODON_LENGTH > len(sequence):
        raise ExtractionError(
            "insufficient_downstream_sequence", "Sequence ends before the full start codon"
        )
    if sequence[cds_start : cds_start + START_CODON_LENGTH] != "ATG":
        observed = sequence[cds_start : cds_start + START_CODON_LENGTH]
        raise ExtractionError("non_atg_start", f"Annotated CDS begins with {observed!r}, not 'ATG'")
    if cds_start + START_CODON_LENGTH + DOWNSTREAM_LENGTH > len(sequence):
        available = len(sequence) - (cds_start + START_CODON_LENGTH)
        raise ExtractionError(
            "insufficient_downstream_sequence",
            f"Only {available} of {DOWNSTREAM_LENGTH} required downstream nucleotides are available",
        )

    window_start = cds_start - UPSTREAM_LENGTH
    window_end = cds_start + START_CODON_LENGTH + DOWNSTREAM_LENGTH
    window = sequence[window_start:window_end]
    if len(window) != WINDOW_LENGTH:
        raise ExtractionError(
            "unexpected_sequence_length",
            f"Expected {WINDOW_LENGTH} nt, extracted {len(window)} nt",
        )
    if any(base not in VALID_BASES for base in window):
        invalid = sorted(set(window) - VALID_BASES)
        raise ExtractionError("invalid_sequence", f"Window contains invalid nucleotide symbols: {invalid}")
    if window[UPSTREAM_LENGTH : UPSTREAM_LENGTH + START_CODON_LENGTH] != "ATG":
        raise ExtractionError("non_atg_start", "Window positions 24–26 do not contain ATG")
    if len(window[:UPSTREAM_LENGTH]) != UPSTREAM_LENGTH:
        raise ExtractionError("insufficient_upstream_sequence", "Upstream region is not 23 nt")
    if len(window[UPSTREAM_LENGTH + START_CODON_LENGTH :]) != DOWNSTREAM_LENGTH:
        raise ExtractionError("insufficient_downstream_sequence", "Downstream region is not 24 nt")

    gc_percent = 100.0 * (window.count("G") + window.count("C")) / WINDOW_LENGTH
    return transcript, gene_name, window, gc_percent


def process_records(path: Path) -> tuple[list[dict[str, str]], list[dict[str, str]], int]:
    """Process a compressed GenBank collection into successes and errors."""
    if not path.is_file():
        raise FileNotFoundError(f"GenBank input not found: {path}")

    successes: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    total_records = 0
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            total_records += 1
            cds_features = [feature for feature in record.features if feature.type == "CDS"]
            transcript = transcript_accession(record)
            gene_name = gene_name_for_record(record, cds_features)
            try:
                transcript, gene_name, sequence, gc_percent = extract_window(record, cds_features)
                successes.append(
                    {
                        "gene_name": gene_name,
                        "transcript": transcript,
                        "sequence": sequence,
                        "GC": f"{gc_percent:.1f}",
                    }
                )
            except ExtractionError as error:
                errors.append(
                    {
                        "transcript": transcript,
                        "gene_name": gene_name,
                        "error_type": error.error_type,
                        "description": str(error),
                    }
                )
                LOG.warning("%s (%s): %s", transcript, error.error_type, error)
    return successes, errors, total_records


def atomic_write(path: Path, writer_function: Callable[[TextIO], None]) -> None:
    """Write one output file through a temporary file before replacing it."""
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as handle:
        temporary_path = Path(handle.name)
        writer_function(handle)
        handle.flush()
        os.fsync(handle.fileno())
    temporary_path.replace(path)


def write_fasta(rows: list[dict[str, str]], handle: TextIO) -> None:
    """Write one FASTA entry per successful transcript, keyed by accession."""
    records = [
        SeqRecord(Seq(row["sequence"]), id=row["transcript"], description="")
        for row in rows
    ]
    SeqIO.write(records, handle, "fasta")


def write_tsv(path: Path, rows: list[dict[str, str]], columns: tuple[str, ...]) -> None:
    """Write tabular results with a stable header and atomic replacement."""
    def writer(handle: TextIO) -> None:
        output = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        output.writeheader()
        output.writerows(rows)

    atomic_write(path, writer)


def gc_distribution(rows: list[dict[str, str]]) -> None:
    """Log summary statistics and a compact histogram of GC percentages."""
    values = [float(row["GC"]) for row in rows]
    if not values:
        LOG.info("GC content distribution: no successful sequences")
        return
    LOG.info(
        "GC%% distribution: min=%.1f, median=%.1f, mean=%.1f, max=%.1f",
        min(values),
        median(values),
        mean(values),
        max(values),
    )
    bins = ((0, 40), (40, 50), (50, 60), (60, 101))
    for lower, upper in bins:
        count = sum(lower <= value < upper for value in values)
        label = f"{lower}–<{upper}%" if upper != 101 else f"{lower}%+"
        LOG.info("GC%% %s: %d", label, count)


def report_validation(
    successes: list[dict[str, str]], errors: list[dict[str, str]], total_records: int
) -> None:
    """Report extraction counts, failures, ATG validation, GC, and reference check."""
    error_counts = Counter(row["error_type"] for row in errors)
    atg_count = sum(row["sequence"][23:26] == "ATG" for row in successes)
    LOG.info("Total GenBank records processed: %d", total_records)
    LOG.info("Successfully extracted sequences: %d", len(successes))
    LOG.info("Failed transcripts: %d", len(errors))
    for error_type, count in sorted(error_counts.items()):
        LOG.info("Failures %s: %d", error_type, count)
    LOG.info("Sequences with valid ATG at positions 24–26: %d", atg_count)
    gc_distribution(successes)

    reference = next(
        (row for row in successes if row["transcript"] == REFERENCE_TRANSCRIPT), None
    )
    if reference is None:
        LOG.warning("Reference transcript %s was not successfully extracted", REFERENCE_TRANSCRIPT)
    elif reference["sequence"] == REFERENCE_SEQUENCE and round(float(reference["GC"])) == 54:
        LOG.info("Reference example %s matches the expected 50-nt sequence and 54%% GC", REFERENCE_TRANSCRIPT)
    else:
        LOG.warning(
            "Reference example mismatch for %s: sequence=%s, GC=%s%%",
            REFERENCE_TRANSCRIPT,
            reference["sequence"],
            reference["GC"],
        )


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")
    try:
        DATA_DIR.mkdir(parents=True, exist_ok=True)
        successes, errors, total_records = process_records(INPUT_FILE)
        atomic_write(FASTA_FILE, lambda handle: write_fasta(successes, handle))
        write_tsv(RESULTS_FILE, successes, ("gene_name", "transcript", "sequence", "GC"))
        write_tsv(
            ERRORS_FILE,
            errors,
            ("transcript", "gene_name", "error_type", "description"),
        )
        report_validation(successes, errors, total_records)
        LOG.info("Wrote FASTA: %s", FASTA_FILE)
        LOG.info("Wrote successful results: %s", RESULTS_FILE)
        LOG.info("Wrote extraction errors: %s", ERRORS_FILE)
    except (OSError, EOFError, gzip.BadGzipFile, ValueError) as error:
        LOG.error("Sequence extraction failed: %s", error)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
