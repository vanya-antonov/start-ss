#!/usr/bin/env python3
"""Extract canonical protein-coding transcripts and transcript-relative CDS coordinates."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


OUTPUT_COLUMNS = [
    "gene_id",
    "gene_name",
    "gene_biotype",
    "transcript_id",
    "transcript_biotype",
    "protein_id",
    "transcript_length",
    "cds_start",
    "cds_end",
    "cds_length",
]
PROGRESS_INTERVAL = 100_000


class GTFFormatError(ValueError):
    """Raised when a GTF record cannot be interpreted safely."""


@dataclass(frozen=True)
class GTFRecord:
    """The GTF fields needed by this analysis."""

    chromosome: str
    feature: str
    start: int
    end: int
    strand: str
    attributes: dict
    line_number: int


def parse_attributes(attribute_text, line_number):
    """Parse GTF attributes, preserving repeated keys such as ``tag``."""
    attributes = defaultdict(list)
    position = 0
    text_length = len(attribute_text)

    while position < text_length:
        while position < text_length and attribute_text[position] in " ;\t\r\n":
            position += 1
        if position >= text_length:
            break

        key_start = position
        while position < text_length and attribute_text[position] not in " ;\t\r\n":
            position += 1
        key = attribute_text[key_start:position]
        while position < text_length and attribute_text[position].isspace():
            position += 1
        if not key or position >= text_length:
            raise GTFFormatError(f"Malformed attribute at GTF line {line_number}")

        if attribute_text[position] == '"':
            position += 1
            value_chars = []
            while position < text_length:
                char = attribute_text[position]
                if char == "\\" and position + 1 < text_length:
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
            while position < text_length and attribute_text[position] != ";":
                position += 1
            value = attribute_text[value_start:position].strip()

        while position < text_length and attribute_text[position].isspace():
            position += 1
        if position < text_length:
            if attribute_text[position] != ";":
                raise GTFFormatError(
                    f"Malformed attribute separator at GTF line {line_number}"
                )
            position += 1
        attributes[key].append(value)

    return dict(attributes)


def first_attribute(attributes, key, default=""):
    values = attributes.get(key, [])
    return values[0] if values else default


def versioned_id(attributes, id_key, version_key):
    """Return a versioned Ensembl ID without appending an existing suffix."""
    identifier = first_attribute(attributes, id_key)
    version = first_attribute(attributes, version_key)
    if identifier and version and not identifier.rsplit(".", maxsplit=1)[-1].isdigit():
        identifier = f"{identifier}.{version}"
    return identifier


def iter_gtf_records(gtf_path, pass_name, canonical_count):
    """Yield relevant records from a gzipped GTF, reporting every 100,000 lines."""
    try:
        with gzip.open(gtf_path, "rt", encoding="utf-8") as gtf_file:
            print(f"{pass_name}: reading {gtf_path.name}", flush=True)
            line_number = 0
            for line_number, line in enumerate(gtf_file, start=1):
                if not line.strip() or line.startswith("#"):
                    if line_number % PROGRESS_INTERVAL == 0:
                        print_progress(pass_name, line_number, canonical_count())
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 9:
                    raise GTFFormatError(
                        f"Expected 9 tab-separated columns at GTF line {line_number}; "
                        f"found {len(fields)}"
                    )
                feature = fields[2]
                if feature not in {"gene", "transcript", "exon", "CDS"}:
                    if line_number % PROGRESS_INTERVAL == 0:
                        print_progress(pass_name, line_number, canonical_count())
                    continue
                try:
                    start, end = int(fields[3]), int(fields[4])
                except ValueError as error:
                    raise GTFFormatError(
                        f"Non-integer coordinates at GTF line {line_number}"
                    ) from error
                if start < 1 or end < start:
                    raise GTFFormatError(
                        f"Invalid 1-based inclusive interval {start}-{end} "
                        f"at GTF line {line_number}"
                    )
                attributes = parse_attributes(fields[8], line_number)
                yield GTFRecord(
                    chromosome=fields[0],
                    feature=feature,
                    start=start,
                    end=end,
                    strand=fields[6],
                    attributes=attributes,
                    line_number=line_number,
                )
                if line_number % PROGRESS_INTERVAL == 0:
                    print_progress(pass_name, line_number, canonical_count())
            if line_number % PROGRESS_INTERVAL != 0:
                print_progress(pass_name, line_number, canonical_count())
    except (OSError, EOFError, UnicodeError) as error:
        raise GTFFormatError(f"Could not read gzipped GTF {gtf_path}: {error}") from error


def print_progress(pass_name, line_count, canonical_count):
    print(
        f"{pass_name}: processed {line_count:,} GTF lines; "
        f"{canonical_count:,} canonical transcripts found",
        flush=True,
    )


def get_transcript_metadata(gtf_path):
    """Collect canonical transcripts and select protein-coding entries."""
    genes = {}
    canonical = {}
    for record in iter_gtf_records(gtf_path, "Pass 1/2", lambda: len(canonical)):
        attributes = record.attributes
        if record.feature == "gene":
            gene_id = versioned_id(attributes, "gene_id", "gene_version")
            if gene_id:
                genes[gene_id] = {
                    "gene_name": first_attribute(attributes, "gene_name"),
                    "gene_biotype": first_attribute(attributes, "gene_biotype"),
                }
        elif record.feature == "transcript":
            tags = [tag for value in attributes.get("tag", []) for tag in value.split(",")]
            if "Ensembl_canonical" not in tags:
                continue
            transcript_id = versioned_id(
                attributes, "transcript_id", "transcript_version"
            )
            gene_id = versioned_id(attributes, "gene_id", "gene_version")
            if not transcript_id or not gene_id:
                raise GTFFormatError(
                    f"Canonical transcript record at GTF line {record.line_number} "
                    "is missing gene_id or transcript_id"
                )
            metadata = {
                "gene_id": gene_id,
                "gene_name": first_attribute(attributes, "gene_name"),
                "gene_biotype": first_attribute(attributes, "gene_biotype"),
                "transcript_id": transcript_id,
                "transcript_biotype": first_attribute(
                    attributes, "transcript_biotype"
                ),
                "protein_id": versioned_id(
                    attributes, "protein_id", "protein_version"
                ),
                "chromosome": record.chromosome,
                "strand": record.strand,
            }
            if transcript_id in canonical and canonical[transcript_id] != metadata:
                raise GTFFormatError(
                    f"Conflicting canonical transcript records for {transcript_id}"
                )
            canonical[transcript_id] = metadata

    selected = {}
    for transcript_id, metadata in canonical.items():
        gene_metadata = genes.get(metadata["gene_id"], {})
        gene_biotype = gene_metadata.get("gene_biotype") or metadata["gene_biotype"]
        if gene_biotype != "protein_coding":
            continue
        if metadata["transcript_biotype"] != "protein_coding":
            continue
        metadata["gene_biotype"] = gene_biotype
        metadata["gene_name"] = metadata["gene_name"] or gene_metadata.get(
            "gene_name", ""
        )
        selected[transcript_id] = metadata

    excluded_count = len(canonical) - len(selected)
    print(
        f"Pass 1/2 complete: {len(canonical):,} canonical transcripts found; "
        f"{len(selected):,} protein-coding transcripts selected; "
        f"{excluded_count:,} non-protein-coding transcripts excluded.",
        flush=True,
    )
    return len(canonical), excluded_count, selected


def collect_transcript_segments(gtf_path, selected_transcripts):
    """Collect exon, CDS, and protein ID records for selected transcripts."""
    exons = defaultdict(list)
    cds_segments = defaultdict(list)
    for record in iter_gtf_records(
        gtf_path, "Pass 2/2", lambda: len(selected_transcripts)
    ):
        if record.feature not in {"exon", "CDS"}:
            continue
        transcript_id = versioned_id(
            record.attributes, "transcript_id", "transcript_version"
        )
        if not transcript_id:
            raise GTFFormatError(
                f"{record.feature} record at GTF line {record.line_number} "
                "is missing transcript_id"
            )
        if transcript_id not in selected_transcripts:
            continue
        segment = (record.chromosome, record.strand, record.start, record.end)
        target = exons if record.feature == "exon" else cds_segments
        target[transcript_id].append(segment)
        if record.feature == "CDS":
            protein_id = versioned_id(
                record.attributes, "protein_id", "protein_version"
            )
            if protein_id and not selected_transcripts[transcript_id]["protein_id"]:
                selected_transcripts[transcript_id]["protein_id"] = protein_id
    print(
        f"Pass 2/2 complete: exon/CDS records collected for "
        f"{len(selected_transcripts):,} canonical protein-coding transcripts.",
        flush=True,
    )
    return exons, cds_segments


def map_cds_to_transcript(metadata, exon_segments, cds_segments):
    """Map CDS intervals onto 1-based coordinates in the spliced transcript."""
    if not exon_segments:
        return None, "missing exon records"
    if not cds_segments:
        return None, "missing CDS records"
    chromosome = metadata["chromosome"]
    strand = metadata["strand"]
    if strand not in {"+", "-"}:
        return None, f"unsupported transcript strand {strand!r}"

    all_segments = exon_segments + cds_segments
    if any(
        seq != chromosome or segment_strand != strand
        for seq, segment_strand, _, _ in all_segments
    ):
        return None, "transcript, exon, or CDS chromosome/strand mismatch"

    ordered_exons = sorted(
        exon_segments, key=lambda segment: segment[2], reverse=strand == "-"
    )
    exon_offsets = []
    offset = 0
    genomic_exons = sorted((start, end) for _, _, start, end in exon_segments)
    if any(
        next_start <= previous_end
        for (_, previous_end), (next_start, _) in zip(
            genomic_exons, genomic_exons[1:]
        )
    ):
        return None, "overlapping exon records"
    for _, _, exon_start, exon_end in ordered_exons:
        exon_offsets.append((exon_start, exon_end, offset))
        offset += exon_end - exon_start + 1

    mapped = []
    for _, _, cds_start, cds_end in cds_segments:
        containing_exons = [
            (exon_start, exon_end, exon_offset)
            for exon_start, exon_end, exon_offset in exon_offsets
            if exon_start <= cds_start and cds_end <= exon_end
        ]
        if len(containing_exons) != 1:
            return None, "CDS segment is not contained in exactly one exon"
        exon_start, exon_end, exon_offset = containing_exons[0]
        if strand == "+":
            transcript_start = exon_offset + cds_start - exon_start + 1
            transcript_end = exon_offset + cds_end - exon_start + 1
        else:
            transcript_start = exon_offset + exon_end - cds_end + 1
            transcript_end = exon_offset + exon_end - cds_start + 1
        mapped.append((transcript_start, transcript_end))

    mapped.sort()
    for (_, previous_end), (next_start, _) in zip(mapped, mapped[1:]):
        if next_start <= previous_end:
            return None, "overlapping CDS segments in transcript coordinates"
        if next_start != previous_end + 1:
            return None, "CDS is not contiguous in spliced transcript coordinates"

    cds_start = mapped[0][0]
    cds_end = mapped[-1][1]
    cds_length = sum(end - start + 1 for start, end in mapped)
    transcript_length = sum(end - start + 1 for _, _, start, end in exon_segments)
    if not (1 <= cds_start <= cds_end <= transcript_length):
        return None, "CDS coordinates fall outside the spliced transcript"
    if cds_length != cds_end - cds_start + 1:
        return None, "CDS length does not match its coordinate span"
    return (cds_start, cds_end, cds_length), None


def build_rows(gtf_path):
    """Build one output row for every canonical protein-coding transcript."""
    canonical_count, excluded_count, selected = get_transcript_metadata(gtf_path)
    exons, cds_segments = collect_transcript_segments(gtf_path, selected)
    rows = []
    cds_issues = []

    for transcript_id, metadata in selected.items():
        exon_segments = exons.get(transcript_id, [])
        cds_coordinates, reason = map_cds_to_transcript(
            metadata, exon_segments, cds_segments.get(transcript_id, [])
        )
        transcript_length = (
            sum(end - start + 1 for _, _, start, end in exon_segments)
            if exon_segments
            else "-"
        )
        if reason:
            cds_issues.append((metadata["gene_id"], transcript_id, reason))
            cds_start = cds_end = cds_length = "-"
        else:
            cds_start, cds_end, cds_length = cds_coordinates

        rows.append(
            [
                metadata["gene_id"] or "-",
                metadata["gene_name"] or "-",
                metadata["gene_biotype"] or "-",
                metadata["transcript_id"] or "-",
                metadata["transcript_biotype"] or "-",
                metadata["protein_id"] or "-",
                transcript_length,
                cds_start,
                cds_end,
                cds_length,
            ]
        )

    rows.sort(key=lambda row: (row[0], row[3]))
    cds_issues.sort(key=lambda row: (row[0], row[1]))
    return canonical_count, excluded_count, rows, cds_issues


def output_path_for_gtf(gtf_path):
    """Place a same-basename canonical transcript TSV beside the input GTF."""
    suffix = ".gtf.gz"
    if not gtf_path.name.endswith(suffix):
        raise ValueError(f"Expected a .gtf.gz input filename: {gtf_path.name}")
    output_name = gtf_path.name[: -len(suffix)] + ".canonical_transcripts.tsv"
    return gtf_path.parent / output_name


def write_text_atomically(path, rows):
    """Write the output TSV atomically, using '-' for missing values."""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as output_file:
            writer = csv.writer(output_file, delimiter="\t", lineterminator="\n")
            writer.writerow(OUTPUT_COLUMNS)
            writer.writerows(
                ["-" if value is None or value == "" else value for value in row]
                for row in rows
            )
        Path(temporary_name).replace(path)
    except Exception:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract canonical protein-coding transcripts and transcript-relative "
            "CDS coordinates from a gzipped Ensembl GTF."
        )
    )
    parser.add_argument("gtf", type=Path, help="Path to an Ensembl .gtf.gz file")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        return 2
    args = parser.parse_args()
    if not args.gtf.is_file():
        parser.error(f"GTF file does not exist: {args.gtf}")

    print(f"Input GTF: {args.gtf}", flush=True)
    try:
        canonical_count, excluded_count, rows, cds_issues = build_rows(args.gtf)
        output_path = output_path_for_gtf(args.gtf)
        print(f"Writing output TSV: {output_path}", flush=True)
        write_text_atomically(output_path, rows)
    except (GTFFormatError, OSError, UnicodeError, ValueError) as error:
        parser.error(f"Could not process GTF {args.gtf}: {error}")

    issue_counts = Counter(reason for _, _, reason in cds_issues)
    print(f"Canonical transcripts found: {canonical_count:,}", flush=True)
    print(f"Canonical protein-coding transcripts written: {len(rows):,}", flush=True)
    print(
        f"Canonical transcripts excluded (not protein-coding): {excluded_count:,}",
        flush=True,
    )
    print(
        f"Transcripts retained with missing or invalid CDS coordinates: "
        f"{len(cds_issues):,}",
        flush=True,
    )
    for reason, count in sorted(issue_counts.items()):
        print(f"  {reason}: {count:,}", flush=True)
    print(f"Output: {output_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
