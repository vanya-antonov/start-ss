#!/usr/bin/env python3
"""Add Ensembl gene and canonical transcript IDs from Ensembl GTF files."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


EXPERIMENT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = EXPERIMENT_DIR.parents[2]
DEFAULT_GTF_DIR = PROJECT_ROOT / "data" / "large_files" / "ensembl" / "gtf"
OUTPUT_COLUMNS = (
    "ensembl_gene",
    "ensembl_canonical_transcript",
    "ensembl_protein",
    "cds_start",
    "cds_end",
    "cds_length",
)

sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "mammals"))
from constants import SPECIES_DICT  # noqa: E402


class GTFFormatError(ValueError):
    """Raised when a GTF record cannot be interpreted safely."""


def parse_attributes(attribute_text: str, line_number: int) -> dict[str, list[str]]:
    """Parse GTF attributes, preserving repeated keys such as ``tag``."""
    attributes: dict[str, list[str]] = defaultdict(list)
    position = 0

    while position < len(attribute_text):
        while position < len(attribute_text) and attribute_text[position] in " ;\t\r\n":
            position += 1
        if position >= len(attribute_text):
            break

        key_start = position
        while position < len(attribute_text) and attribute_text[position] not in " ;\t\r\n":
            position += 1
        key = attribute_text[key_start:position]
        while position < len(attribute_text) and attribute_text[position].isspace():
            position += 1
        if not key or position >= len(attribute_text):
            raise GTFFormatError(f"Malformed attribute at GTF line {line_number}")

        if attribute_text[position] == '"':
            position += 1
            value_chars = []
            while position < len(attribute_text):
                char = attribute_text[position]
                if char == "\\" and position + 1 < len(attribute_text):
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
            while position < len(attribute_text) and attribute_text[position] != ";":
                position += 1
            value = attribute_text[value_start:position].strip()

        while position < len(attribute_text) and attribute_text[position].isspace():
            position += 1
        if position < len(attribute_text):
            if attribute_text[position] != ";":
                raise GTFFormatError(
                    f"Malformed attribute separator at GTF line {line_number}"
                )
            position += 1
        attributes[key].append(value)

    return dict(attributes)


def first_attribute(attributes: dict[str, list[str]], key: str) -> str:
    """Return the first value for an attribute, or an empty string."""
    values = attributes.get(key, [])
    return values[0] if values else ""


def versioned_id(
    attributes: dict[str, list[str]], id_key: str, version_key: str
) -> str:
    """Return an Ensembl ID with its version, whether embedded or separate."""
    identifier = first_attribute(attributes, id_key)
    version = first_attribute(attributes, version_key)
    if identifier and version and not identifier.rsplit(".", maxsplit=1)[-1].isdigit():
        identifier = f"{identifier}.{version}"
    return identifier


def gtf_path_for_species(species_name: str, gtf_dir: Path) -> Path:
    """Find the release 116 GTF using the scientific-name filename prefix."""
    if species_name not in SPECIES_DICT:
        raise ValueError(
            f"Species key {species_name!r} is not present in SPECIES_DICT."
        )
    scientific_name = str(SPECIES_DICT[species_name].get("species", "")).strip()
    if not scientific_name:
        raise ValueError(f"Species key {species_name!r} has no scientific name.")
    filename_prefix = scientific_name.replace(" ", "_")
    filename_prefix = filename_prefix[0].upper() + filename_prefix[1:]
    matches = sorted(gtf_dir.glob(f"{filename_prefix}.*.116.gtf.gz"))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected one release 116 GTF for {species_name!r} matching "
            f"{filename_prefix}.*.116.gtf.gz in {gtf_dir}; found {len(matches)}."
        )
    return matches[0]


def read_input_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read and validate the ortholog TSV."""
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


@dataclass(frozen=True)
class CanonicalTranscript:
    """GTF metadata for one canonical transcript."""

    gene_id: str
    transcript_id: str
    gene_name: str
    gene_biotype: str
    transcript_biotype: str
    protein_id: str
    chromosome: str
    strand: str


@dataclass(frozen=True)
class AnnotationResult:
    """Resolved identifiers, CDS coordinates, and status for one gene symbol."""

    gene_id: str = "-"
    transcript_id: str = "-"
    protein_id: str = "-"
    cds_start: str = "-"
    cds_end: str = "-"
    cds_length: str = "-"
    gene_status: str = "missing"
    transcript_status: str = "not_assessed"
    protein_status: str = "not_assessed"
    cds_status: str = "not_assessed"
    reason: str = ""


GTFSegment = tuple[str, str, int, int]


def read_species_gtf(
    gtf_path: Path, requested_symbols: set[str]
) -> tuple[
    dict[str, dict[str, str]],
    dict[str, dict[str, CanonicalTranscript]],
    dict[str, list[GTFSegment]],
    dict[str, list[GTFSegment]],
    dict[str, set[str]],
]:
    """Collect matching genes, canonical transcripts, exons, and CDS in one pass."""
    genes_by_symbol: dict[str, dict[str, str]] = defaultdict(dict)
    canonical_by_gene: dict[str, dict[str, CanonicalTranscript]] = defaultdict(dict)
    exons_by_transcript: dict[str, list[GTFSegment]] = defaultdict(list)
    cds_by_transcript: dict[str, list[GTFSegment]] = defaultdict(list)
    proteins_by_transcript: dict[str, set[str]] = defaultdict(set)
    matching_gene_ids: set[str] = set()
    line_count = 0

    try:
        with gzip.open(gtf_path, "rt", encoding="utf-8") as gtf_file:
            for line_count, line in enumerate(gtf_file, start=1):
                if line_count % 100_000 == 0:
                    report_gtf_progress(
                        gtf_path,
                        line_count,
                        len(matching_gene_ids),
                        sum(len(items) for items in canonical_by_gene.values()),
                    )
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 9:
                    raise GTFFormatError(
                        f"Expected 9 tab-separated columns at GTF line {line_count}; "
                        f"found {len(fields)}"
                    )
                feature = fields[2]
                if feature not in {"gene", "transcript", "exon", "CDS"}:
                    continue

                attributes = parse_attributes(fields[8], line_count)
                gene_name = first_attribute(attributes, "gene_name")
                gene_id = versioned_id(attributes, "gene_id", "gene_version")
                if feature == "gene":
                    if gene_name in requested_symbols and gene_id:
                        genes_by_symbol[gene_name][gene_id] = first_attribute(
                            attributes, "gene_biotype"
                        )
                        matching_gene_ids.add(gene_id)
                elif feature == "transcript":
                    if gene_name not in requested_symbols and gene_id not in matching_gene_ids:
                        continue
                    tags = {
                        tag
                        for value in attributes.get("tag", [])
                        for tag in value.split(",")
                    }
                    if "Ensembl_canonical" not in tags:
                        continue
                    transcript_id = versioned_id(
                        attributes, "transcript_id", "transcript_version"
                    )
                    if not gene_id or not transcript_id:
                        raise GTFFormatError(
                            f"Canonical transcript at GTF line {line_count} is "
                            "missing gene_id or transcript_id"
                        )
                    metadata = CanonicalTranscript(
                        gene_id=gene_id,
                        transcript_id=transcript_id,
                        gene_name=gene_name,
                        gene_biotype=first_attribute(attributes, "gene_biotype"),
                        transcript_biotype=first_attribute(
                            attributes, "transcript_biotype"
                        ),
                        protein_id=versioned_id(
                            attributes, "protein_id", "protein_version"
                        ),
                        chromosome=fields[0],
                        strand=fields[6],
                    )
                    existing = canonical_by_gene[gene_id].get(transcript_id)
                    if existing and existing != metadata:
                        raise GTFFormatError(
                            f"Conflicting canonical transcript records for "
                            f"{transcript_id}"
                        )
                    canonical_by_gene[gene_id][transcript_id] = metadata
                    if metadata.protein_id:
                        proteins_by_transcript[transcript_id].add(metadata.protein_id)
                elif feature in {"exon", "CDS"}:
                    if gene_name not in requested_symbols and gene_id not in matching_gene_ids:
                        continue
                    transcript_id = versioned_id(
                        attributes, "transcript_id", "transcript_version"
                    )
                    if not transcript_id:
                        raise GTFFormatError(
                            f"{feature} record at GTF line {line_count} "
                            "is missing transcript_id"
                        )
                    try:
                        start, end = int(fields[3]), int(fields[4])
                    except ValueError as error:
                        raise GTFFormatError(
                            f"Non-integer coordinates at GTF line {line_count}"
                        ) from error
                    if start < 1 or end < start:
                        raise GTFFormatError(
                            f"Invalid 1-based inclusive interval {start}-{end} "
                            f"at GTF line {line_count}"
                        )
                    segment = (fields[0], fields[6], start, end)
                    if feature == "exon":
                        exons_by_transcript[transcript_id].append(segment)
                    else:
                        cds_by_transcript[transcript_id].append(segment)
                        protein_id = versioned_id(
                            attributes, "protein_id", "protein_version"
                        )
                        if protein_id:
                            proteins_by_transcript[transcript_id].add(protein_id)

    except (OSError, EOFError, UnicodeError) as error:
        raise GTFFormatError(f"Could not read gzipped GTF {gtf_path}: {error}") from error

    report_gtf_progress(
        gtf_path,
        line_count,
        len(matching_gene_ids),
        sum(len(items) for items in canonical_by_gene.values()),
        complete=True,
    )
    return (
        genes_by_symbol,
        canonical_by_gene,
        exons_by_transcript,
        cds_by_transcript,
        proteins_by_transcript,
    )


def report_gtf_progress(
    gtf_path: Path,
    line_count: int,
    gene_count: int,
    transcript_count: int,
    complete: bool = False,
) -> None:
    """Print periodic or completion progress for one GTF scan."""
    stage = "GTF scan complete" if complete else "GTF progress"
    print(
        f"{stage} ({gtf_path.name}): {line_count:,} lines; "
        f"{gene_count:,} matching gene(s); "
        f"{transcript_count:,} canonical transcript(s).",
        flush=True,
    )


def map_cds_to_transcript(
    transcript: CanonicalTranscript,
    exons: list[GTFSegment],
    cds_segments: list[GTFSegment],
) -> tuple[tuple[int, int, int] | None, str, str]:
    """Map CDS genomic intervals to 1-based inclusive spliced-transcript coordinates."""
    if not exons:
        return None, "missing", "missing exon records"
    if not cds_segments:
        return None, "missing", "missing CDS records"
    if transcript.strand not in {"+", "-"}:
        return None, "invalid", f"unsupported transcript strand {transcript.strand!r}"

    all_segments = exons + cds_segments
    if any(
        chromosome != transcript.chromosome or strand != transcript.strand
        for chromosome, strand, _, _ in all_segments
    ):
        return None, "invalid", "transcript, exon, or CDS chromosome/strand mismatch"

    genomic_exons = sorted((start, end) for _, _, start, end in exons)
    if any(
        next_start <= previous_end
        for (_, previous_end), (next_start, _) in zip(
            genomic_exons, genomic_exons[1:]
        )
    ):
        return None, "invalid", "overlapping exon records"

    ordered_exons = sorted(
        exons,
        key=lambda segment: segment[2],
        reverse=transcript.strand == "-",
    )
    exon_offsets = []
    transcript_length = 0
    for _, _, exon_start, exon_end in ordered_exons:
        exon_offsets.append((exon_start, exon_end, transcript_length))
        transcript_length += exon_end - exon_start + 1

    mapped_cds = []
    for _, _, cds_start, cds_end in cds_segments:
        containing_exons = [
            (exon_start, exon_end, offset)
            for exon_start, exon_end, offset in exon_offsets
            if exon_start <= cds_start and cds_end <= exon_end
        ]
        if len(containing_exons) != 1:
            return None, "invalid", "CDS segment is not contained in exactly one exon"
        exon_start, exon_end, offset = containing_exons[0]
        if transcript.strand == "+":
            tx_start = offset + cds_start - exon_start + 1
            tx_end = offset + cds_end - exon_start + 1
        else:
            tx_start = offset + exon_end - cds_end + 1
            tx_end = offset + exon_end - cds_start + 1
        mapped_cds.append((tx_start, tx_end))

    mapped_cds.sort()
    for (_, previous_end), (next_start, _) in zip(mapped_cds, mapped_cds[1:]):
        if next_start <= previous_end:
            return None, "invalid", "overlapping CDS segments in transcript coordinates"
        if next_start != previous_end + 1:
            return None, "invalid", "CDS is not contiguous in spliced transcript coordinates"

    cds_start = mapped_cds[0][0]
    cds_end = mapped_cds[-1][1]
    cds_length = sum(end - start + 1 for start, end in mapped_cds)
    if not (1 <= cds_start <= cds_end <= transcript_length):
        return None, "invalid", "CDS coordinates fall outside the spliced transcript"
    if cds_length != cds_end - cds_start + 1:
        return None, "invalid", "CDS length does not match its coordinate span"
    return (cds_start, cds_end, cds_length), "valid", ""


def resolve_symbol(
    symbol: str,
    genes_by_symbol: dict[str, dict[str, str]],
    canonical_by_gene: dict[str, dict[str, CanonicalTranscript]],
    exons_by_transcript: dict[str, list[GTFSegment]],
    cds_by_transcript: dict[str, list[GTFSegment]],
    proteins_by_transcript: dict[str, set[str]],
) -> AnnotationResult:
    """Resolve one exact gene symbol without choosing among ambiguous matches."""
    gene_matches = genes_by_symbol.get(symbol, {})
    if not symbol or symbol == "-" or not gene_matches:
        return AnnotationResult(reason="no gene with an exact gene_name match")
    if len(gene_matches) > 1:
        return AnnotationResult(
            gene_status="ambiguous",
            reason=f"multiple gene IDs match gene_name {symbol!r}",
        )

    gene_id, gene_biotype = next(iter(gene_matches.items()))
    candidates = {}
    for transcript_id, metadata in canonical_by_gene.get(gene_id, {}).items():
        selected_gene_biotype = gene_biotype or metadata.gene_biotype
        if (
            selected_gene_biotype == "protein_coding"
            and metadata.transcript_biotype == "protein_coding"
        ):
            candidates[transcript_id] = metadata
    if not candidates:
        return AnnotationResult(
            gene_id=gene_id,
            gene_status="matched",
            transcript_status="missing",
            reason="no canonical protein-coding transcript",
        )
    if len(candidates) > 1:
        return AnnotationResult(
            gene_id=gene_id,
            gene_status="matched",
            transcript_status="ambiguous",
            reason="multiple canonical protein-coding transcripts match the gene",
        )

    transcript_id, metadata = next(iter(candidates.items()))
    protein_ids = set(proteins_by_transcript.get(transcript_id, set()))
    if metadata.protein_id:
        protein_ids.add(metadata.protein_id)
    if len(protein_ids) == 1:
        protein_id = next(iter(protein_ids))
        protein_status = "matched"
        protein_reason = ""
    elif not protein_ids:
        protein_id = "-"
        protein_status = "missing"
        protein_reason = "no protein_id in canonical transcript or CDS records"
    else:
        protein_id = "-"
        protein_status = "ambiguous"
        protein_reason = "multiple distinct protein IDs are associated with the transcript"

    coordinates, cds_status, cds_reason = map_cds_to_transcript(
        metadata,
        exons_by_transcript.get(transcript_id, []),
        cds_by_transcript.get(transcript_id, []),
    )
    if coordinates is None:
        cds_start = cds_end = cds_length = "-"
    else:
        cds_start, cds_end, cds_length = (str(value) for value in coordinates)

    reasons = [reason for reason in (protein_reason, cds_reason) if reason]
    return AnnotationResult(
        gene_id=gene_id,
        transcript_id=transcript_id,
        protein_id=protein_id,
        cds_start=cds_start,
        cds_end=cds_end,
        cds_length=cds_length,
        gene_status="matched",
        transcript_status="matched",
        protein_status=protein_status,
        cds_status=cds_status,
        reason="; ".join(reasons),
    )


def write_output(
    path: Path, input_columns: list[str], rows: list[dict[str, str]]
) -> None:
    """Write all input rows and append the Ensembl result columns."""
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
            {
                column: row.get(column)
                if row.get(column) not in {None, ""}
                else "-"
                for column in columns
            }
            for row in rows
        )


def automatic_paths(gene_name: str) -> tuple[Path, Path]:
    """Build input and output paths from the required gene name."""
    safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", gene_name.strip()).strip("._-")
    if not safe_name:
        raise ValueError("Gene name must contain filename-safe characters")
    return (
        EXPERIMENT_DIR / f"01_{safe_name}_ncbi_orthologs.tsv",
        EXPERIMENT_DIR / f"02_{safe_name}_ncbi_orthologs.tsv",
    )


def process_table(input_path: Path, output_path: Path, gtf_dir: Path) -> int:
    """Resolve gene, canonical transcript, protein, and CDS data for each row."""
    print(f"Reading input table: {input_path}", flush=True)
    input_columns, rows = read_input_table(input_path)
    print(f"Input rows: {len(rows)}", flush=True)
    rows_by_species: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        species_name = row["species_name"]
        if species_name not in SPECIES_DICT:
            raise ValueError(
                f"Species key {species_name!r} is not present in SPECIES_DICT."
            )
        rows_by_species[species_name].append(row)

    overall = Counter()
    for species_name, species_rows in rows_by_species.items():
        gtf_path = gtf_path_for_species(species_name, gtf_dir)
        symbols = {
            row.get("ortholog_gene_symbol") or ""
            for row in species_rows
            if (row.get("ortholog_gene_symbol") or "") not in {"", "-"}
        }
        print(
            f"\nProcessing species {species_name} ({len(species_rows)} row(s), "
            f"{len(symbols)} unique symbol(s)).",
            flush=True,
        )
        print(f"Reading GTF: {gtf_path.name}", flush=True)
        (
            genes_by_symbol,
            canonical_by_gene,
            exons_by_transcript,
            cds_by_transcript,
            proteins_by_transcript,
        ) = read_species_gtf(gtf_path, symbols)

        matched_gene_count = sum(
            len(genes_by_symbol.get(symbol, {})) for symbol in symbols
        )
        canonical_count = sum(
            1
            for transcripts in canonical_by_gene.values()
            for metadata in transcripts.values()
            if metadata.transcript_biotype == "protein_coding"
        )
        print(
            f"Matched gene records: {matched_gene_count}; "
            f"canonical protein-coding transcripts: {canonical_count}.",
            flush=True,
        )

        results = {
            symbol: resolve_symbol(
                symbol,
                genes_by_symbol,
                canonical_by_gene,
                exons_by_transcript,
                cds_by_transcript,
                proteins_by_transcript,
            )
            for symbol in symbols
        }
        species_counts = Counter()
        for row_number, row in enumerate(species_rows, start=1):
            symbol = row.get("ortholog_gene_symbol") or ""
            result = results.get(symbol) or resolve_symbol(
                symbol,
                genes_by_symbol,
                canonical_by_gene,
                exons_by_transcript,
                cds_by_transcript,
                proteins_by_transcript,
            )
            row.update(
                {
                    "ensembl_gene": result.gene_id,
                    "ensembl_canonical_transcript": result.transcript_id,
                    "ensembl_protein": result.protein_id,
                    "cds_start": result.cds_start,
                    "cds_end": result.cds_end,
                    "cds_length": result.cds_length,
                }
            )
            species_counts["rows"] += 1
            species_counts[f"gene_{result.gene_status}"] += 1
            species_counts[f"transcript_{result.transcript_status}"] += 1
            species_counts[f"protein_{result.protein_status}"] += 1
            species_counts[f"cds_{result.cds_status}"] += 1
            if (
                result.gene_status == "matched"
                and result.transcript_status == "matched"
                and result.protein_status == "matched"
                and result.cds_status == "valid"
            ):
                species_counts["complete_matches"] += 1
            if result.reason:
                print(
                    f"{species_name} row {row_number} ({symbol!r}): "
                    f"{result.reason}.",
                    flush=True,
                )

        print(f"Species summary ({species_name}):", flush=True)
        print(f"  Rows: {species_counts['rows']}", flush=True)
        print(f"  Complete matches: {species_counts['complete_matches']}", flush=True)
        print(f"  Missing genes: {species_counts['gene_missing']}", flush=True)
        print(f"  Ambiguous gene matches: {species_counts['gene_ambiguous']}", flush=True)
        print(
            f"  Missing canonical transcripts: {species_counts['transcript_missing']}",
            flush=True,
        )
        print(
            f"  Ambiguous canonical transcripts: {species_counts['transcript_ambiguous']}",
            flush=True,
        )
        print(f"  Missing CDS annotations: {species_counts['cds_missing']}", flush=True)
        print(f"  Invalid CDS annotations: {species_counts['cds_invalid']}", flush=True)
        print(f"  Missing protein IDs: {species_counts['protein_missing']}", flush=True)
        print(f"  Ambiguous protein IDs: {species_counts['protein_ambiguous']}", flush=True)
        overall.update(species_counts)

    print("\nWriting output table.", flush=True)
    write_output(output_path, input_columns, rows)
    print("\nFinal summary across species:", flush=True)
    print(f"  Input rows: {overall['rows']}", flush=True)
    print(f"  Complete matches: {overall['complete_matches']}", flush=True)
    print(f"  Missing genes: {overall['gene_missing']}", flush=True)
    print(f"  Ambiguous gene matches: {overall['gene_ambiguous']}", flush=True)
    print(
        f"  Missing canonical transcripts: {overall['transcript_missing']}",
        flush=True,
    )
    print(
        f"  Ambiguous canonical transcripts: {overall['transcript_ambiguous']}",
        flush=True,
    )
    print(f"  Missing CDS annotations: {overall['cds_missing']}", flush=True)
    print(f"  Invalid CDS annotations: {overall['cds_invalid']}", flush=True)
    print(f"  Missing protein IDs: {overall['protein_missing']}", flush=True)
    print(f"  Ambiguous protein IDs: {overall['protein_ambiguous']}", flush=True)
    print(f"  Output: {output_path}", flush=True)
    return 0


def main() -> int:
    """Parse paths and add Ensembl IDs to the GABARAP ortholog table."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "gene_name",
        help="Gene name used to select the experiment input/output filenames",
    )
    parser.add_argument("--input", type=Path, help="Override the input TSV path")
    parser.add_argument("--output", type=Path, help="Override the output TSV path")
    parser.add_argument(
        "--gtf-dir",
        type=Path,
        default=DEFAULT_GTF_DIR,
        help="Directory with Ensembl 116 GTF files",
    )
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        return 2
    args = parser.parse_args()
    try:
        default_input, default_output = automatic_paths(args.gene_name)
        input_path = args.input or default_input
        output_path = args.output or default_output
        return process_table(input_path, output_path, args.gtf_dir)
    except (OSError, UnicodeError, csv.Error, GTFFormatError, ValueError) as error:
        parser.error(f"Could not add canonical transcript IDs: {error}")


if __name__ == "__main__":
    raise SystemExit(main())
