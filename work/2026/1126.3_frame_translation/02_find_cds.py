#!/usr/bin/env python3
"""Find the first two ATG-initiated CDS candidates in each FASTA record."""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


STOP_CODONS = {"TAA", "TAG", "TGA"}
OUTPUT_COLUMNS = (
    "name",
    "cds_1_start",
    "cds_1_end",
    "cds_1_length",
    "cds_2_start",
    "cds_2_end",
    "cds_2_length",
)
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


@dataclass
class CdsCandidate:
    """A start codon and its first in-frame stop codon, if present."""

    start: int | None
    end: int | None
    length: int | None
    nucleotide_sequence: str | None


@dataclass
class TranscriptCds:
    """CDS candidate results associated with a complete FASTA header."""

    name: str
    sequence: str
    cds_1: CdsCandidate
    cds_2: CdsCandidate


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Find the first two ATG-started CDS candidates in each transcript, "
            "using each start codon's own reading frame."
        )
    )
    parser.add_argument("--fasta", required=True, type=Path, help="Input multi-record FASTA file")
    parser.add_argument("--output", required=True, type=Path, help="Output CDS coordinates TSV")
    parser.add_argument(
        "--cds_2_faa",
        type=Path,
        help="Optional FASTA output for complete CDS 2 protein translations",
    )
    return parser.parse_args()


def read_fasta_records(fasta_path: Path) -> Iterator[tuple[str, str]]:
    """Yield full FASTA headers and uppercase sequences from a multi-record file."""
    if not fasta_path.is_file():
        raise FileNotFoundError(f"Input FASTA file does not exist: {fasta_path}")

    header: str | None = None
    sequence_parts: list[str] = []
    line_count = 0
    record_count = 0
    print(f"Reading multi-record FASTA: {fasta_path}")
    with fasta_path.open(encoding="utf-8") as fasta_file:
        for line_count, line in enumerate(fasta_file, start=1):
            line = line.rstrip("\r\n")
            if not line.strip():
                continue
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(sequence_parts).upper()
                    if not sequence:
                        raise ValueError(f"FASTA record {header!r} has an empty sequence.")
                    record_count += 1
                    yield header, sequence
                header = line[1:]
                if not header.strip():
                    raise ValueError(f"FASTA header is empty at line {line_count}.")
                sequence_parts = []
            elif header is None:
                raise ValueError(f"Sequence data found before the first FASTA header at line {line_count}.")
            else:
                sequence_parts.append("".join(line.split()))

            if line_count % 100_000 == 0:
                print(f"  Read {line_count:,} FASTA lines; {record_count:,} records completed")

    if header is None:
        raise ValueError(f"No FASTA records found in {fasta_path}.")
    sequence = "".join(sequence_parts).upper()
    if not sequence:
        raise ValueError(f"FASTA record {header!r} has an empty sequence.")
    record_count += 1
    yield header, sequence
    print(f"  Loaded {record_count:,} FASTA records")


def find_first_two_start_codons(sequence: str) -> tuple[int | None, int | None]:
    """Return the 1-based positions of the first and second ATG occurrences."""
    first_index = sequence.find("ATG")
    if first_index == -1:
        return None, None
    second_index = sequence.find("ATG", first_index + 1)
    return first_index + 1, (second_index + 1 if second_index != -1 else None)


def find_in_frame_stop(sequence: str, start_position: int) -> int | None:
    """Return the 1-based inclusive end of the first in-frame stop codon."""
    start_index = start_position - 1
    for codon_index in range(start_index + 3, len(sequence) - 2, 3):
        if sequence[codon_index : codon_index + 3] in STOP_CODONS:
            return codon_index + 3
    return None


def find_candidate(sequence: str, start_position: int | None) -> CdsCandidate:
    """Build one CDS candidate from a start position, retaining incomplete starts."""
    if start_position is None:
        return CdsCandidate(None, None, None, None)

    end_position = find_in_frame_stop(sequence, start_position)
    if end_position is None:
        return CdsCandidate(start_position, None, None, None)

    length = end_position - start_position + 1
    nucleotide_sequence = sequence[start_position - 1 : end_position]
    return CdsCandidate(start_position, end_position, length, nucleotide_sequence)


def translate_cds(nucleotide_sequence: str) -> str:
    """Translate a complete CDS and omit its terminal stop symbol."""
    amino_acids = "".join(
        CODON_TABLE.get(nucleotide_sequence[index : index + 3], "X")
        for index in range(0, len(nucleotide_sequence), 3)
    )
    if amino_acids.endswith("*"):
        amino_acids = amino_acids[:-1]
    return amino_acids


def format_value(value: int | None) -> int | str:
    return value if value is not None else "-"


def write_tsv(results: list[TranscriptCds], output_path: Path) -> None:
    """Write CDS coordinates and lengths with the requested column order."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing CDS coordinate TSV: {output_path}")
    with output_path.open("w", newline="", encoding="utf-8") as tsv_file:
        writer = csv.writer(tsv_file, delimiter="\t", lineterminator="\n")
        writer.writerow(OUTPUT_COLUMNS)
        for result in results:
            writer.writerow((
                result.name,
                format_value(result.cds_1.start),
                format_value(result.cds_1.end),
                format_value(result.cds_1.length),
                format_value(result.cds_2.start),
                format_value(result.cds_2.end),
                format_value(result.cds_2.length),
            ))


def write_protein_fasta(results: list[TranscriptCds], output_path: Path) -> None:
    """Write translations only for complete CDS 2 candidates."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    complete_candidates = [
        result for result in results if result.cds_2.nucleotide_sequence is not None
    ]
    print(
        f"Writing {len(complete_candidates):,} complete CDS 2 protein translations: "
        f"{output_path}"
    )
    with output_path.open("w", encoding="utf-8") as protein_file:
        for result in complete_candidates:
            protein = translate_cds(result.cds_2.nucleotide_sequence or "")
            protein_file.write(f">{result.name}\n")
            for line_start in range(0, len(protein), 60):
                protein_file.write(f"{protein[line_start : line_start + 60]}\n")


def validate_output_paths(fasta_path: Path, output_path: Path, protein_path: Path | None) -> None:
    """Prevent output paths from overwriting an input or each other."""
    input_resolved = fasta_path.resolve()
    output_resolved = output_path.resolve()
    if output_resolved == input_resolved:
        raise ValueError("TSV output path must differ from the input FASTA path.")
    if protein_path is not None:
        protein_resolved = protein_path.resolve()
        if protein_resolved in {input_resolved, output_resolved}:
            raise ValueError("Protein FASTA output path must differ from input and TSV output paths.")


def main() -> int:
    arguments = parse_args()
    try:
        validate_output_paths(arguments.fasta, arguments.output, arguments.cds_2_faa)
        results: list[TranscriptCds] = []
        total_records = 0
        complete_cds_1 = 0
        complete_cds_2 = 0

        for record_number, (header, sequence) in enumerate(read_fasta_records(arguments.fasta), start=1):
            total_records += 1
            print(f"Processing transcript {record_number:,}: {header}")
            first_start, second_start = find_first_two_start_codons(sequence)
            cds_1 = find_candidate(sequence, first_start)
            cds_2 = find_candidate(sequence, second_start)
            complete_cds_1 += cds_1.nucleotide_sequence is not None
            complete_cds_2 += cds_2.nucleotide_sequence is not None
            print(
                f"  CDS 1: start={format_value(cds_1.start)}, "
                f"end={format_value(cds_1.end)}, length={format_value(cds_1.length)}"
            )
            print(
                f"  CDS 2: start={format_value(cds_2.start)}, "
                f"end={format_value(cds_2.end)}, length={format_value(cds_2.length)}"
            )
            results.append(TranscriptCds(header, sequence, cds_1, cds_2))

        write_tsv(results, arguments.output)
        if arguments.cds_2_faa is not None:
            write_protein_fasta(results, arguments.cds_2_faa)

        print(f"Processed {total_records:,} transcripts")
        print(f"  Complete CDS 1 candidates: {complete_cds_1:,}")
        print(f"  Complete CDS 2 candidates: {complete_cds_2:,}")
        print(f"  TSV output: {arguments.output}")
        if arguments.cds_2_faa is not None:
            print(f"  CDS 2 protein FASTA: {arguments.cds_2_faa}")
    except (OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
