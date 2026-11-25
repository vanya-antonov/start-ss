#!/usr/bin/env python3
"""Create an interactive Plotly view of three translated transcript frames."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import plotly.graph_objects as go


SCRIPT_DIR = Path(__file__).resolve().parent
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
    species_name: str
    transcript_id: str
    sequence: str
    cds_start: int | None
    cds_end: int | None


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """Read a multiline FASTA while reporting progress for large inputs."""
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence_parts: list[str] = []
    line_count = 0

    print(f"Reading transcript FASTA: {path}")
    with path.open(encoding="utf-8") as fasta_file:
        for line_count, line in enumerate(fasta_file, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence_parts)))
                header = line[1:].split(maxsplit=1)[0]
                sequence_parts = []
            elif header is None:
                raise ValueError(f"Sequence data found before the first FASTA header at line {line_count}.")
            else:
                sequence_parts.append(line)
            if line_count % 100_000 == 0:
                print(f"  Read {line_count:,} FASTA lines; {len(records):,} records completed")

    if header is not None:
        records.append((header, "".join(sequence_parts)))
    if not records:
        raise ValueError(f"No FASTA records found in {path}.")
    print(f"  Loaded {len(records):,} transcript records")
    return records


def read_annotations(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    """Index ortholog rows by species and transcript ID."""
    annotations: dict[tuple[str, str], dict[str, str]] = {}
    print(f"Reading CDS annotations: {path}")
    with path.open(newline="", encoding="utf-8") as tsv_file:
        reader = csv.DictReader(tsv_file, delimiter="\t")
        required_columns = {"species_name", "ensembl_canonical_transcript", "cds_start", "cds_end"}
        missing = required_columns.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Annotation TSV is missing required columns: {', '.join(sorted(missing))}")
        for row in reader:
            key = (row["species_name"].strip(), row["ensembl_canonical_transcript"].strip())
            if all(key):
                annotations[key] = row
    print(f"  Loaded {len(annotations):,} transcript annotation rows")
    return annotations


def parse_header(header: str) -> tuple[str, str]:
    """Parse the expected species_name|transcript_id FASTA identifier."""
    if "|" in header:
        species_name, transcript_id = header.split("|", 1)
        return species_name, transcript_id
    return "unknown", header


def parse_coordinate(value: str | None) -> int | None:
    if value is None or value.strip() in {"", "-"}:
        return None
    try:
        coordinate = int(value)
    except ValueError:
        return None
    return coordinate if coordinate > 0 else None


def make_transcripts(
    fasta_records: list[tuple[str, str]],
    annotations: dict[tuple[str, str], dict[str, str]],
) -> list[Transcript]:
    transcripts: list[Transcript] = []
    annotation_matches = 0
    fallback_count = 0
    for header, sequence in fasta_records:
        species_name, transcript_id = parse_header(header)
        row = annotations.get((species_name, transcript_id))
        cds_start = parse_coordinate(row.get("cds_start") if row else None)
        cds_end = parse_coordinate(row.get("cds_end") if row else None)
        if cds_start is not None:
            annotation_matches += 1
        else:
            fallback_count += 1
        transcripts.append(
            Transcript(species_name, transcript_id, sequence.upper().replace("U", "T"), cds_start, cds_end)
        )
    print(f"  Matched {annotation_matches:,} transcripts to a CDS start; {fallback_count:,} use transcript base 1 as fallback")
    return transcripts


def iter_codons(sequence: str, offset: int) -> Iterator[tuple[int, str, str]]:
    """Yield zero-based codon starts, codons, and one-letter amino acids."""
    for start in range(offset, len(sequence) - 2, 3):
        codon = sequence[start : start + 3]
        yield start, codon, CODON_TABLE.get(codon, "X")


def add_transcript_traces(fig: go.Figure, transcript: Transcript, transcript_index: int) -> list[float]:
    """Add three frame tracks and their codon/marker overlays to the figure."""
    y_values: list[float] = []
    base_y = transcript_index * 4.0
    sequence_length = len(transcript.sequence)

    # Keep all transcript nucleotides visible. Frame 0 is the phase whose codon
    # boundaries include cds_start; frames +1/+2 shift that phase by one/two nt.
    reference_offset = ((transcript.cds_start - 1) % 3) if transcript.cds_start is not None and transcript.cds_start <= sequence_length else 0
    if transcript.cds_start is not None and transcript.cds_start > sequence_length:
        print(f"  Warning: CDS start {transcript.cds_start} is beyond {transcript.species_name}|{transcript.transcript_id}; using transcript base 1")
    for frame in range(3):
        y = base_y + frame
        y_values.append(y)
        offset = (reference_offset + frame) % 3
        codons = list(iter_codons(transcript.sequence, offset))
        if not codons:
            continue

        line_x: list[float | None] = []
        line_y: list[float | None] = []
        aa_x: list[float] = []
        aa_y: list[float] = []
        amino_acids: list[str] = []
        hover_text: list[str] = []
        start_x: list[float] = []
        start_hover: list[str] = []
        stop_x: list[float] = []
        stop_hover: list[str] = []

        for codon_start, codon, amino_acid in codons:
            first_base = codon_start + 1  # transcript coordinates are displayed as 1-based
            last_base = codon_start + 3
            midpoint = codon_start + 2
            line_x.extend((first_base, last_base, None))
            line_y.extend((y, y, None))
            aa_x.append(midpoint)
            aa_y.append(y)
            amino_acids.append(amino_acid)
            hover_text.append(f"{transcript.species_name}|{transcript.transcript_id}<br>Frame {frame}<br>nt {first_base}–{last_base}: {codon}<br>Amino acid: {amino_acid}")
            if codon in START_CODONS:
                start_x.append(midpoint)
                start_hover.append(f"{transcript.species_name}|{transcript.transcript_id}<br>Frame {frame}<br>Start codon {codon} at nt {first_base}–{last_base}")
            if codon in STOP_CODONS:
                stop_x.append(midpoint)
                stop_hover.append(f"{transcript.species_name}|{transcript.transcript_id}<br>Frame {frame}<br>Stop codon {codon} at nt {first_base}–{last_base}")

        color = FRAME_COLORS[frame]
        frame_label = ("frame 0", "frame +1", "frame +2")[frame]
        fig.add_trace(go.Scatter(
            x=line_x, y=line_y, mode="lines", name=frame_label,
            legendgroup=frame_label, showlegend=transcript_index == 0,
            line=dict(color=color, width=10), hoverinfo="skip",
        ))
        fig.add_trace(go.Scatter(
            x=aa_x, y=aa_y, mode="text", text=amino_acids,
            textposition="middle center", textfont=dict(color="white", size=9),
            customdata=hover_text, hovertemplate="%{customdata}<extra></extra>",
            name=f"{frame_label} amino acids", legendgroup=frame_label, showlegend=False,
        ))
        if start_x:
            fig.add_trace(go.Scatter(
                x=start_x, y=[y] * len(start_x), mode="markers", text=start_hover,
                marker=dict(symbol="line-ew", size=14, color="white", line=dict(color="#333333", width=1)),
                name="Start codon (ATG)", legendgroup="start", showlegend=False,
                hovertemplate="%{text}<extra></extra>",
            ))
        if stop_x:
            fig.add_trace(go.Scatter(
                x=stop_x, y=[y] * len(stop_x), mode="markers", text=stop_hover,
                marker=dict(symbol="line-ew", size=14, color="black", line=dict(color="white", width=1)),
                name="Stop codon (TAA/TAG/TGA)", legendgroup="stop", showlegend=False,
                hovertemplate="%{text}<extra></extra>",
            ))

    # Shade the annotated CDS interval behind frame 0. Coordinates in the TSV
    # are 1-based inclusive, matching the displayed transcript nucleotide axis.
    if (
        transcript.cds_start is not None
        and transcript.cds_end is not None
        and transcript.cds_start <= transcript.cds_end <= sequence_length
    ):
        fig.add_shape(
            type="rect", xref="x", yref="y",
            x0=transcript.cds_start - 0.5, x1=transcript.cds_end + 0.5,
            y0=base_y - 0.25, y1=base_y + 0.25,
            fillcolor="rgba(229, 57, 53, 0.22)", line=dict(color="rgba(229, 57, 53, 0.7)", width=1),
            layer="below",
        )
    elif transcript.cds_start is not None or transcript.cds_end is not None:
        print(f"  Warning: invalid/incomplete CDS interval for {transcript.species_name}|{transcript.transcript_id}; CDS shading omitted")
    return y_values


def create_figure(transcripts: list[Transcript]) -> go.Figure:
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=[None], y=[None], mode="markers", name="Start codon (ATG)",
        marker=dict(symbol="line-ew", size=14, color="white", line=dict(color="#333333", width=1)),
        hoverinfo="skip", legendgroup="start",
    ))
    fig.add_trace(go.Scatter(
        x=[None], y=[None], mode="markers", name="Stop codon (TAA/TAG/TGA)",
        marker=dict(symbol="line-ew", size=14, color="black", line=dict(color="white", width=1)),
        hoverinfo="skip", legendgroup="stop",
    ))
    tick_values: list[float] = []
    tick_labels: list[str] = []
    max_length = 1
    for index, transcript in enumerate(transcripts):
        print(f"  Plotting {index + 1:,}/{len(transcripts):,}: {transcript.species_name}|{transcript.transcript_id}")
        y_values = add_transcript_traces(fig, transcript, index)
        tick_values.extend(y_values)
        tick_labels.extend(
            f"{transcript.species_name}|{transcript.transcript_id} · {label}"
            for label in ("frame 0", "frame +1", "frame +2")
        )
        max_length = max(max_length, len(transcript.sequence))

    fig.update_layout(
        title="Three-frame transcript translations",
        xaxis=dict(title="Transcript nucleotide position (1-based)", range=[0.5, max_length + 0.5], showgrid=True, zeroline=False),
        yaxis=dict(title="Species | transcript ID · reading frame", tickmode="array", tickvals=tick_values, ticktext=tick_labels, autorange="reversed", showgrid=False, zeroline=False),
        height=max(600, 62 * len(transcripts) + 130),
        hovermode="closest",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
        margin=dict(l=320, r=35, t=115, b=70),
        plot_bgcolor="white",
        paper_bgcolor="white",
    )
    fig.add_annotation(
        xref="paper", yref="paper", x=1, y=-0.14, xanchor="right", showarrow=False,
        text="Frame 0 phase is anchored at annotated cds_start; the full transcript is shown. CDS interval is shaded. Missing coordinates use transcript nucleotide 1 as the frame reference.",
        font=dict(size=10, color="#555555"),
    )
    return fig


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Translate every transcript in three reading frames and create an interactive Plotly HTML page.",
        epilog="Example: python work/2026/1125.3_frame_plot/01_plot_transcript_frames.py work/2026/1125.3_frame_plot/00_GABARAP_transcripts.fasta",
    )
    parser.add_argument("fasta", type=Path, help="Input transcript FASTA; headers should be species_name|transcript_id")
    parser.add_argument(
        "--annotations", type=Path, default=SCRIPT_DIR / "00_GABARAP_orthologs.tsv",
        help="TSV with species_name, ensembl_canonical_transcript, cds_start, and cds_end (default: %(default)s)",
    )
    parser.add_argument(
        "--output", type=Path, default=SCRIPT_DIR / "01_GABARAP_transcript_frames.html",
        help="Output standalone HTML file (default: %(default)s)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.fasta.is_file():
        raise SystemExit(f"Input FASTA not found: {args.fasta}")
    if not args.annotations.is_file():
        raise SystemExit(f"Annotation TSV not found: {args.annotations}")

    fasta_records = read_fasta(args.fasta)
    annotations = read_annotations(args.annotations)
    transcripts = make_transcripts(fasta_records, annotations)
    print("Building Plotly figure")
    figure = create_figure(transcripts)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.write_html(args.output, include_plotlyjs=True, full_html=True)
    print(f"Wrote interactive transcript plot for {len(transcripts):,} transcripts to {args.output}")


if __name__ == "__main__":
    main()
