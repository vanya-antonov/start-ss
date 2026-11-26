#!/usr/bin/env python3
"""Plot transcript-level CDS annotations and translations in three frames."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import plotly.graph_objects as go


ANNOTATION_DIR = Path("data/large_files/ensembl/gtf")
FRAME_COLORS = ("#e53935", "#2e9d50", "#2878d0")
START_CODONS = {"ATG"}
STOP_CODONS = {"TAA", "TAG", "TGA"}
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
class Transcript:
    label: str
    transcript_id: str
    sequence: str
    cds_start: int | None
    cds_end: int | None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create an interactive Plotly HTML page with an annotated CDS track "
            "and translations in three transcript-based reading frames."
        )
    )
    parser.add_argument("fasta", type=Path, help="Input transcript FASTA file")
    parser.add_argument(
        "--annotations-dir",
        type=Path,
        default=ANNOTATION_DIR,
        help=f"Directory containing release 116 canonical transcript TSVs (default: {ANNOTATION_DIR})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output HTML path (default: 01_<species>_3_frames.html beside the FASTA)",
    )
    return parser.parse_args()


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield FASTA identifiers and sequences, reporting large-file progress."""
    header: str | None = None
    sequence_parts: list[str] = []
    line_count = 0
    record_count = 0
    print(f"Reading transcript FASTA: {path}")
    with path.open(encoding="utf-8") as fasta_file:
        for line_count, line in enumerate(fasta_file, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    record_count += 1
                    yield header, "".join(sequence_parts)
                header = line[1:].split(maxsplit=1)[0]
                sequence_parts = []
            elif header is None:
                raise ValueError(f"Sequence data found before the first FASTA header at line {line_count}.")
            else:
                sequence_parts.append(line)
            if line_count % 100_000 == 0:
                print(f"  Read {line_count:,} lines; {record_count:,} transcript records completed")
    if header is not None:
        record_count += 1
        yield header, "".join(sequence_parts)
    print(f"  Loaded {record_count:,} transcript records")


def read_annotations(directory: Path) -> dict[str, tuple[int | None, int | None]]:
    """Load transcript-relative CDS coordinates from release 116 TSV files."""
    annotation_files = sorted(directory.glob("*.116.canonical_transcripts.tsv"))
    if not annotation_files:
        raise FileNotFoundError(f"No *.116.canonical_transcripts.tsv files found in {directory}")

    annotations: dict[str, tuple[int | None, int | None]] = {}
    print(f"Reading CDS coordinates from {len(annotation_files):,} annotation TSV files in {directory}")
    for path in annotation_files:
        print(f"  Reading {path}")
        with path.open(newline="", encoding="utf-8") as tsv_file:
            reader = csv.DictReader(tsv_file, delimiter="\t")
            required = {"transcript_id", "cds_start", "cds_end"}
            missing = required.difference(reader.fieldnames or [])
            if missing:
                raise ValueError(f"{path} is missing required columns: {', '.join(sorted(missing))}")
            row_count = 0
            for row_count, row in enumerate(reader, start=1):
                transcript_id = row["transcript_id"].strip()
                if transcript_id:
                    annotations[transcript_id] = (
                        parse_coordinate(row.get("cds_start")),
                        parse_coordinate(row.get("cds_end")),
                    )
                if row_count % 100_000 == 0:
                    print(f"    Read {row_count:,} annotation rows from {path.name}")
    print(f"  Indexed {len(annotations):,} transcript annotations")
    return annotations


def parse_coordinate(value: str | None) -> int | None:
    if value is None or value.strip() in {"", "-"}:
        return None
    try:
        coordinate = int(value)
    except ValueError:
        return None
    return coordinate if coordinate > 0 else None


def parse_fasta_id(header: str) -> tuple[str, str]:
    """Return a display label and the transcript ID from species|transcript IDs."""
    if "|" in header:
        species_name, transcript_id = header.split("|", 1)
        return f"{species_name}|{transcript_id}", transcript_id
    return header, header


def make_transcripts(
    fasta_path: Path,
    annotations: dict[str, tuple[int | None, int | None]],
) -> list[Transcript]:
    transcripts: list[Transcript] = []
    matched = 0
    for header, sequence in iter_fasta(fasta_path):
        label, transcript_id = parse_fasta_id(header)
        cds_start, cds_end = annotations.get(transcript_id, (None, None))
        if cds_start is not None and cds_end is not None:
            matched += 1
        transcripts.append(
            Transcript(label, transcript_id, sequence.upper().replace("U", "T"), cds_start, cds_end)
        )
    print(f"  Matched CDS coordinates for {matched:,} of {len(transcripts):,} transcripts")
    return transcripts


def add_transcript_traces(
    figure: go.Figure,
    transcript: Transcript,
    transcript_index: int,
    first_transcript: bool,
) -> None:
    """Add a CDS track and three colored translation tracks for one transcript."""
    sequence_length = len(transcript.sequence)
    valid_cds = (
        transcript.cds_start is not None
        and transcript.cds_end is not None
        and transcript.cds_start <= transcript.cds_end <= sequence_length
    )
    visibility = transcript_index == 0
    if valid_cds:
        cds_start = transcript.cds_start
        cds_end = transcript.cds_end
        assert cds_start is not None and cds_end is not None
        figure.add_trace(go.Scatter(
            x=[cds_start, cds_end], y=[3, 3], mode="lines+markers",
            line=dict(color="#6b7280", width=14),
            marker=dict(color="#374151", size=7),
            name="Annotated CDS", legendgroup="cds", showlegend=first_transcript,
            visible=visibility,
            customdata=[[transcript.label, cds_start, cds_end], [transcript.label, cds_start, cds_end]],
            hovertemplate=("%{customdata[0]}<br>Annotated CDS: nt %{customdata[1]}–%{customdata[2]}"
                           "<extra></extra>"),
        ))
    else:
        print(f"  Warning: CDS coordinates missing or invalid for {transcript.label}; CDS track omitted")

    for frame in range(3):
        codons = [
            (start, transcript.sequence[start : start + 3])
            for start in range(frame, sequence_length - 2, 3)
        ]
        x_values = [start + 2 for start, _ in codons]
        amino_acids = [CODON_TABLE.get(codon, "X") for _, codon in codons]
        hover = [
            f"{transcript.label}<br>Frame {frame:+d}<br>nt {start + 1}–{start + 3}: {codon}"
            f"<br>Amino acid: {CODON_TABLE.get(codon, 'X')}"
            for start, codon in codons
        ]
        figure.add_trace(go.Scatter(
            x=x_values, y=[2 - frame] * len(codons), mode="markers+text",
            marker=dict(color=FRAME_COLORS[frame], size=18, symbol="square"),
            text=amino_acids, textposition="middle center",
            textfont=dict(color="white", size=9), customdata=hover,
            hovertemplate="%{customdata}<extra></extra>",
            name=f"Frame {frame}", legendgroup=f"frame{frame}",
            showlegend=first_transcript, visible=visibility,
        ))

        start_codons = [(start + 2, codon) for start, codon in codons if codon in START_CODONS]
        stop_codons = [(start + 2, codon) for start, codon in codons if codon in STOP_CODONS]
        for codon_type, matches, color, symbol in (
            ("Start", start_codons, "#ffffff", "triangle-up"),
            ("Stop", stop_codons, "#111827", "x"),
        ):
            if not matches:
                continue
            figure.add_trace(go.Scatter(
                x=[position for position, _ in matches], y=[2 - frame] * len(matches),
                mode="markers",
                marker=dict(color=color, symbol=symbol, size=13, line=dict(color="#111827", width=1)),
                customdata=[
                    f"{transcript.label}<br>Frame {frame}<br>{codon_type} codon {codon} at nt {position - 1}–{position + 1}"
                    for position, codon in matches
                ],
                hovertemplate="%{customdata}<extra></extra>",
                name=f"{codon_type} codon", legendgroup=codon_type.lower(),
                showlegend=first_transcript and frame == 0, visible=visibility,
            ))


def build_figure(transcripts: list[Transcript]) -> go.Figure:
    figure = go.Figure()
    trace_ranges: list[tuple[int, int]] = []
    for transcript_index, transcript in enumerate(transcripts):
        print(f"Preparing plot {transcript_index + 1:,}/{len(transcripts):,}: {transcript.label}")
        first_trace = len(figure.data)
        add_transcript_traces(figure, transcript, transcript_index, transcript_index == 0)
        trace_ranges.append((first_trace, len(figure.data)))

    buttons = []
    for selected_index, transcript in enumerate(transcripts):
        first_trace, last_trace = trace_ranges[selected_index]
        visible = [
            first_trace <= trace_index < last_trace
            for trace_index in range(len(figure.data))
        ]
        buttons.append({
            "label": transcript.label,
            "method": "update",
            "args": [
                {"visible": visible},
                {"title": f"{transcript.label} — CDS and three-frame translation",
                 "xaxis.range": [0, len(transcript.sequence) + 1]},
            ],
        })

    first = transcripts[0]
    figure.update_layout(
        title=f"{first.label} — CDS and three-frame translation",
        xaxis=dict(title="Transcript nucleotide position (1-based)", range=[0, len(first.sequence) + 1],
                   showgrid=True, dtick=100),
        yaxis=dict(
            tickmode="array", tickvals=[3, 2, 1, 0],
            ticktext=["Annotated CDS", "Frame 0", "Frame +1", "Frame +2"],
            range=[-0.6, 3.6], fixedrange=True,
        ),
        updatemenus=[dict(
            buttons=buttons, direction="down", showactive=True,
            x=1, xanchor="right", y=1.14, yanchor="top",
        )],
        height=480,
        hovermode="closest",
        margin=dict(l=110, r=30, t=120, b=70),
        legend=dict(orientation="h", y=-0.2),
    )
    return figure


def main() -> None:
    args = parse_args()
    if not args.fasta.is_file():
        raise FileNotFoundError(f"Input FASTA does not exist: {args.fasta}")
    annotations = read_annotations(args.annotations_dir)
    transcripts = make_transcripts(args.fasta, annotations)
    if not transcripts:
        raise ValueError(f"No FASTA records found in {args.fasta}")
    fasta_stem = args.fasta.stem.removeprefix("00_").removesuffix("_transcripts")
    output = args.output or args.fasta.with_name(f"01_{fasta_stem}_3_frames.html")
    output.parent.mkdir(parents=True, exist_ok=True)
    figure = build_figure(transcripts)
    print(f"Writing interactive Plotly HTML: {output}")
    figure.write_html(output, include_plotlyjs=True, full_html=True)
    print(f"Wrote plot for {len(transcripts):,} transcripts to {output}")


if __name__ == "__main__":
    main()
