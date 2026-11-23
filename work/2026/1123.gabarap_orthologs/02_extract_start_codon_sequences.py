#!/usr/bin/env python3
"""Extract validated 50-nt windows around annotated transcript start codons."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import TextIO


EXPERIMENT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = EXPERIMENT_DIR.parents[2]
DEFAULT_CDNA_DIR = PROJECT_ROOT / "data" / "large_files" / "ensembl" / "cdna"
UPSTREAM_LENGTH = 23
START_CODON_LENGTH = 3
DOWNSTREAM_LENGTH = 24
WINDOW_LENGTH = UPSTREAM_LENGTH + START_CODON_LENGTH + DOWNSTREAM_LENGTH
VALID_BASES = frozenset("ACGT")
PROGRESS_INTERVAL = 100_000
RESULT_COLUMNS = (
    "species_name",
    "species",
    "gene_symbol",
    "ensembl_gene",
    "transcript",
    "ensembl_protein",
    "cds_start",
    "cds_end",
    "cds_length",
    "sequence",
    "GC",
)
ERROR_COLUMNS = (
    "species_name",
    "species",
    "gene_symbol",
    "ensembl_gene",
    "transcript",
    "error_type",
    "description",
)

sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "mammals"))
from constants import SPECIES_DICT  # noqa: E402


class ExtractionError(ValueError):
    """A per-row extraction failure that can be recorded in the errors TSV."""

    def __init__(self, error_type: str, description: str) -> None:
        super().__init__(description)
        self.error_type = error_type


def read_input_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read and validate the annotated ortholog TSV."""
    if not path.is_file():
        raise FileNotFoundError(f"Input TSV not found: {path}")
    with path.open(encoding="utf-8-sig", newline="") as input_file:
        reader = csv.DictReader(input_file, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Input TSV is empty or has no header")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("Input TSV contains duplicate column names")
        required = {
            "species_name",
            "species",
            "ortholog_gene_symbol",
            "ensembl_gene",
            "ensembl_canonical_transcript",
            "ensembl_protein",
            "cds_start",
            "cds_end",
            "cds_length",
        }
        missing = required - set(reader.fieldnames)
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
            if row["species_name"] not in SPECIES_DICT:
                raise ValueError(
                    f"Unknown species key {row['species_name']!r} at line {line_number}"
                )
            rows.append(row)
    return list(reader.fieldnames), rows


def fasta_for_species(species_name: str, cdna_dir: Path) -> Path:
    """Find the release 116 cDNA FASTA from the SPECIES_DICT scientific name."""
    scientific_name = str(SPECIES_DICT[species_name].get("species", "")).strip()
    if not scientific_name:
        raise ValueError(f"Species key {species_name!r} has no scientific name")
    prefix = scientific_name.replace(" ", "_")
    prefix = prefix[0].upper() + prefix[1:]
    matches = sorted(cdna_dir.glob(f"{prefix}.*.cdna.all.fa.gz"))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected one cDNA FASTA for {species_name!r}, matching "
            f"{prefix}.*.cdna.all.fa.gz in {cdna_dir}; found {len(matches)}"
        )
    return matches[0]


def find_requested_sequences(
    fasta_path: Path, requested_ids: set[str]
) -> tuple[dict[str, str], set[str]]:
    """Read one gzipped FASTA and retain only requested versioned transcript IDs."""
    found: dict[str, str] = {}
    duplicates: set[str] = set()
    current_id = ""
    sequence_parts: list[str] = []
    line_count = 0

    def save_current_record() -> None:
        if current_id not in requested_ids:
            return
        sequence = "".join(sequence_parts).upper()
        if current_id in found:
            duplicates.add(current_id)
            found.pop(current_id, None)
        elif current_id not in duplicates:
            found[current_id] = sequence

    print(f"Reading cDNA FASTA: {fasta_path.name}", flush=True)
    try:
        with gzip.open(fasta_path, "rt", encoding="utf-8") as fasta_file:
            for line_count, line in enumerate(fasta_file, start=1):
                if line_count % PROGRESS_INTERVAL == 0:
                    print(
                        f"  Processed {line_count:,} FASTA lines; "
                        f"{len(found):,} requested transcripts found.",
                        flush=True,
                    )
                if line.startswith(">"):
                    save_current_record()
                    header = line[1:].strip()
                    current_id = header.split(maxsplit=1)[0] if header else ""
                    sequence_parts = []
                elif current_id:
                    sequence_parts.append("".join(line.split()))
                elif line.strip():
                    raise ValueError(
                        f"Sequence data before the first FASTA header at line {line_count}"
                    )
            save_current_record()
    except (OSError, EOFError, UnicodeError) as error:
        raise ValueError(f"Could not read gzipped FASTA {fasta_path}: {error}") from error

    print(
        f"  Read {line_count:,} FASTA lines; {len(found):,} requested transcript(s) "
        f"found; {len(duplicates):,} duplicate transcript ID(s).",
        flush=True,
    )
    return found, duplicates


def validate_annotation(row: dict[str, str]) -> tuple[str, int, int, int]:
    """Return required IDs and validated 1-based CDS coordinates for one row."""
    transcript_id = row["ensembl_canonical_transcript"]
    if not transcript_id or transcript_id == "-":
        raise ExtractionError("missing_transcript", "No canonical transcript ID")
    try:
        cds_start, cds_end, cds_length = (
            int(row[column]) for column in ("cds_start", "cds_end", "cds_length")
        )
    except (TypeError, ValueError) as error:
        raise ExtractionError("missing_cds_coordinates", "CDS coordinates are missing or non-numeric") from error
    if cds_start < 1 or cds_end < cds_start or cds_length != cds_end - cds_start + 1:
        raise ExtractionError("invalid_cds_coordinates", "CDS coordinate range or length is invalid")
    return transcript_id, cds_start, cds_end, cds_length


def extract_window(sequence: str, cds_start: int) -> tuple[str, float]:
    """Extract 23 nt + annotated ATG + 24 nt in transcript orientation."""
    start_index = cds_start - 1
    if start_index < UPSTREAM_LENGTH:
        raise ExtractionError(
            "insufficient_upstream_sequence",
            f"CDS starts at transcript position {cds_start}; need {UPSTREAM_LENGTH} upstream nt",
        )
    if start_index + START_CODON_LENGTH > len(sequence):
        raise ExtractionError("insufficient_downstream_sequence", "Transcript ends before the full start codon")
    observed_start = sequence[start_index : start_index + START_CODON_LENGTH]
    if observed_start != "ATG":
        raise ExtractionError(
            "non_atg_start", f"Annotated CDS begins with {observed_start!r}, not 'ATG'"
        )
    window_end = start_index + START_CODON_LENGTH + DOWNSTREAM_LENGTH
    if window_end > len(sequence):
        available = len(sequence) - (start_index + START_CODON_LENGTH)
        raise ExtractionError(
            "insufficient_downstream_sequence",
            f"Only {available} of {DOWNSTREAM_LENGTH} downstream nt are available",
        )
    window = sequence[start_index - UPSTREAM_LENGTH : window_end]
    if len(window) != WINDOW_LENGTH:
        raise ExtractionError(
            "unexpected_sequence_length",
            f"Expected {WINDOW_LENGTH} nt, extracted {len(window)} nt",
        )
    invalid_bases = sorted(set(window) - VALID_BASES)
    if invalid_bases:
        raise ExtractionError(
            "invalid_sequence",
            f"Window contains invalid nucleotide symbols: {invalid_bases}",
        )
    gc_percent = 100.0 * (window.count("G") + window.count("C")) / WINDOW_LENGTH
    return window, gc_percent


def atomic_write(path: Path, writer_function) -> None:
    """Write a result file through a temporary file before replacing it."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as output_file:
        temporary_path = Path(output_file.name)
        writer_function(output_file)
        output_file.flush()
        os.fsync(output_file.fileno())
    temporary_path.replace(path)


def write_tsv(path: Path, rows: list[dict[str, str]], columns: tuple[str, ...]) -> None:
    """Write a TSV with a stable header."""
    def writer(handle: TextIO) -> None:
        output = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        output.writeheader()
        output.writerows(rows)

    atomic_write(path, writer)


def write_fasta(path: Path, rows: list[dict[str, str]]) -> None:
    """Write successful windows as FASTA records keyed by species and transcript."""
    def writer(handle: TextIO) -> None:
        for row in rows:
            record_id = f"{row['species_name']}|{row['transcript']}"
            handle.write(f">{record_id}\n{row['sequence']}\n")

    atomic_write(path, writer)


def automatic_paths(gene_name: str) -> tuple[Path, Path, Path, Path]:
    """Build default input and output paths from the required gene name."""
    safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", gene_name.strip()).strip("._-")
    if not safe_name:
        raise ValueError("Gene name must contain filename-safe characters")
    return (
        EXPERIMENT_DIR / f"01_{safe_name}_orthologs.tsv",
        EXPERIMENT_DIR / f"02_{safe_name}_start_codon_windows.fasta",
        EXPERIMENT_DIR / f"02_{safe_name}_start_codon_windows.tsv",
        EXPERIMENT_DIR / f"02_{safe_name}_extraction_errors.tsv",
    )


def process_table(
    input_path: Path,
    cdna_dir: Path,
    fasta_output: Path,
    results_output: Path,
    errors_output: Path,
) -> int:
    """Extract windows for all valid annotated transcript rows."""
    print(f"Reading annotated ortholog table: {input_path}", flush=True)
    _, rows = read_input_table(input_path)
    print(f"Input rows: {len(rows):,}", flush=True)

    grouped_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped_rows[row["species_name"]].append(row)

    sequences_by_species: dict[str, dict[str, str]] = {}
    duplicates_by_species: dict[str, set[str]] = {}
    annotation_errors: dict[int, ExtractionError] = {}
    requested_by_species: dict[str, set[str]] = defaultdict(set)
    for row_index, row in enumerate(rows):
        try:
            transcript_id, _, _, _ = validate_annotation(row)
            requested_by_species[row["species_name"]].add(transcript_id)
        except ExtractionError as error:
            annotation_errors[row_index] = error

    for species_name, species_rows in grouped_rows.items():
        requested_ids = requested_by_species.get(species_name, set())
        if not requested_ids:
            print(
                f"\nSkipping cDNA FASTA scan for {species_name}: "
                "no rows have valid transcript/CDS annotations.",
                flush=True,
            )
            sequences_by_species[species_name] = {}
            duplicates_by_species[species_name] = set()
            continue
        fasta_path = fasta_for_species(species_name, cdna_dir)
        print(
            f"\nProcessing {species_name}: {len(species_rows)} row(s), "
            f"{len(requested_ids)} transcript ID(s).",
            flush=True,
        )
        sequences, duplicates = find_requested_sequences(fasta_path, requested_ids)
        sequences_by_species[species_name] = sequences
        duplicates_by_species[species_name] = duplicates

    successes: list[dict[str, str]] = []
    errors: list[dict[str, str]] = []
    counts = Counter()
    for row_index, row in enumerate(rows):
        species_name = row["species_name"]
        species = row["species"]
        symbol = row["ortholog_gene_symbol"]
        gene_id = row["ensembl_gene"]
        transcript_id = row["ensembl_canonical_transcript"]
        try:
            if row_index in annotation_errors:
                raise annotation_errors[row_index]
            transcript_id, cds_start, cds_end, cds_length = validate_annotation(row)
            if transcript_id in duplicates_by_species.get(species_name, set()):
                raise ExtractionError(
                    "duplicate_fasta_transcript",
                    f"Transcript ID {transcript_id} occurs more than once in the cDNA FASTA",
                )
            sequence = sequences_by_species.get(species_name, {}).get(transcript_id)
            if sequence is None:
                raise ExtractionError(
                    "transcript_not_found",
                    f"Transcript ID {transcript_id} was not found in the cDNA FASTA",
                )
            if cds_end > len(sequence):
                raise ExtractionError(
                    "cds_out_of_bounds",
                    f"CDS end {cds_end} exceeds transcript length {len(sequence)}",
                )
            window, gc_percent = extract_window(sequence, cds_start)
        except ExtractionError as error:
            errors.append(
                {
                    "species_name": species_name,
                    "species": species,
                    "gene_symbol": symbol,
                    "ensembl_gene": gene_id or "-",
                    "transcript": transcript_id or "-",
                    "error_type": error.error_type,
                    "description": str(error),
                }
            )
            counts[error.error_type] += 1
            print(
                f"{species_name} {symbol}: {error.error_type}: {error}", flush=True
            )
            continue

        successes.append(
            {
                "species_name": species_name,
                "species": species,
                "gene_symbol": symbol,
                "ensembl_gene": gene_id or "-",
                "transcript": transcript_id,
                "ensembl_protein": row["ensembl_protein"] or "-",
                "cds_start": str(cds_start),
                "cds_end": str(cds_end),
                "cds_length": str(cds_length),
                "sequence": window,
                "GC": f"{gc_percent:.1f}",
            }
        )
        counts["success"] += 1

    write_fasta(fasta_output, successes)
    write_tsv(results_output, successes, RESULT_COLUMNS)
    write_tsv(errors_output, errors, ERROR_COLUMNS)
    print(f"\nInput rows: {len(rows):,}", flush=True)
    print(f"Successfully extracted windows: {counts['success']:,}", flush=True)
    print(f"Failed rows: {len(errors):,}", flush=True)
    for error_type, count in sorted(counts.items()):
        if error_type != "success":
            print(f"  {error_type}: {count:,}", flush=True)
    print(f"FASTA: {fasta_output}", flush=True)
    print(f"Results TSV: {results_output}", flush=True)
    print(f"Errors TSV: {errors_output}", flush=True)
    return 0


def main() -> int:
    """Parse arguments and extract transcript start-codon windows."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "gene_name",
        help="Gene name used to select default input and output filenames",
    )
    parser.add_argument("--input", type=Path, help="Override the input TSV path")
    parser.add_argument(
        "--cdna-dir",
        type=Path,
        default=DEFAULT_CDNA_DIR,
        help="Folder containing Ensembl release 116 cDNA FASTA files",
    )
    parser.add_argument("--fasta-output", type=Path, help="Override output FASTA path")
    parser.add_argument("--results-output", type=Path, help="Override successful results TSV path")
    parser.add_argument("--errors-output", type=Path, help="Override errors TSV path")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        return 2
    args = parser.parse_args()
    try:
        input_path, fasta_output, results_output, errors_output = automatic_paths(
            args.gene_name
        )
        return process_table(
            args.input or input_path,
            args.cdna_dir,
            args.fasta_output or fasta_output,
            args.results_output or results_output,
            args.errors_output or errors_output,
        )
    except (OSError, UnicodeError, csv.Error, ValueError) as error:
        parser.error(f"Sequence extraction failed: {error}")


if __name__ == "__main__":
    raise SystemExit(main())
