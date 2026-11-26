#!/usr/bin/env python3
"""Plot one transcript, its annotated CDS, and start/stop codons in three frames."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from html import escape as html_escape
from pathlib import Path
import sys
from typing import Sequence

import plotly.graph_objects as go


NUCLEOTIDE_COLORS = {
    "A": "#F8766D",
    "C": "#00BA38",
    "G": "#619CFF",
    "T": "#619CFF",
}
FRAME_COLORS = ("#F8766D", "#00BA38", "#619CFF")
FRAME_NAMES = ("Frame 0", "Frame +1", "Frame +2")
FRAME_BOUNDS = ((2, 3), (1, 2), (0, 1))
START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a standalone Plotly HTML visualization of one transcript, "
            "its annotated CDS, and start/stop codons in three reading frames."
        )
    )
    parser.add_argument("--fasta", required=True, type=Path, help="FASTA file containing one transcript")
    parser.add_argument("--cds_start", type=int, help="Annotated CDS start (1-based, inclusive)")
    parser.add_argument("--cds_end", type=int, help="Annotated CDS end (1-based, inclusive)")
    parser.add_argument("--output", required=True, type=Path, help="Output standalone HTML file")
    return parser.parse_args()


def read_fasta(fasta_path: Path) -> tuple[str, str]:
    """Read exactly one FASTA record and return its header and sequence."""
    if not fasta_path.is_file():
        raise FileNotFoundError(f"FASTA file does not exist: {fasta_path}")

    transcript_id: str | None = None
    sequence_parts: list[str] = []
    record_count = 0
    line_count = 0
    print(f"Reading FASTA: {fasta_path}")
    with fasta_path.open(encoding="utf-8") as fasta_file:
        for line_count, line in enumerate(fasta_file, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                record_count += 1
                if record_count > 1:
                    raise ValueError("Input FASTA must contain exactly one transcript record.")
                transcript_id = line[1:].strip()
                if not transcript_id:
                    raise ValueError(f"FASTA header is empty at line {line_count}.")
            elif transcript_id is None:
                raise ValueError(f"Sequence data found before the FASTA header at line {line_count}.")
            else:
                sequence_parts.append("".join(line.split()))
            if line_count % 100_000 == 0:
                print(f"  Read {line_count:,} FASTA lines")

    if transcript_id is None:
        raise ValueError(f"No FASTA record found in {fasta_path}.")
    sequence = "".join(sequence_parts).upper()
    if not sequence:
        raise ValueError(f"FASTA record {transcript_id!r} has an empty sequence.")
    print(f"  Loaded {len(sequence):,} nucleotides for {transcript_id}")
    return transcript_id, sequence


def validate_input(
    sequence: str,
    cds_start: int,
    cds_end: int,
) -> tuple[int, int]:
    """Validate 1-based inclusive CDS coordinates and return length and frame."""
    sequence_length = len(sequence)
    if cds_start < 1 or cds_end < 1:
        raise ValueError("CDS coordinates must be positive 1-based positions.")
    if cds_start > cds_end:
        raise ValueError("CDS start must be less than or equal to CDS end.")
    if cds_end > sequence_length:
        raise ValueError(
            f"CDS end ({cds_end}) exceeds transcript length ({sequence_length})."
        )
    cds_length = cds_end - cds_start + 1
    cds_frame = (cds_start - 1) % 3
    return cds_length, cds_frame


def find_cds_coordinates(transcript_id: str) -> tuple[int, int] | None:
    """Search canonical transcript TSVs for this transcript's CDS interval."""
    annotation_dir = Path("data/large_files/ensembl/gtf")
    annotation_files = sorted(annotation_dir.glob("*.canonical_transcripts.tsv"))
    if not annotation_files:
        print(
            f"Error: no *.canonical_transcripts.tsv files found in {annotation_dir}.",
            file=sys.stderr,
        )
        return None

    print(f"Looking up CDS coordinates for {transcript_id} in {len(annotation_files):,} annotation TSV files")
    for annotation_path in annotation_files:
        print(f"  Searching {annotation_path}")
        with annotation_path.open(newline="", encoding="utf-8") as annotation_file:
            reader = csv.DictReader(annotation_file, delimiter="\t")
            required = {"transcript_id", "cds_start", "cds_end"}
            missing_columns = required.difference(reader.fieldnames or [])
            if missing_columns:
                print(
                    f"Error: {annotation_path} is missing columns: {', '.join(sorted(missing_columns))}.",
                    file=sys.stderr,
                )
                continue
            for row_number, row in enumerate(reader, start=1):
                if row_number % 100_000 == 0:
                    print(f"    Read {row_number:,} rows from {annotation_path.name}")
                if row["transcript_id"].strip() != transcript_id:
                    continue
                try:
                    start = int(row["cds_start"])
                    end = int(row["cds_end"])
                except (TypeError, ValueError):
                    return None
                if start <= 0 or end <= 0:
                    return None
                print(f"  Found CDS coordinates {start}–{end} in {annotation_path.name}")
                return start, end
    return None


def identify_codons(sequence: str) -> list[list[tuple[int, str, str]]]:
    """Find start and stop codons independently in each full-transcript frame.

    Each result item is (1-based codon-start position, codon sequence, type).
    """
    codons_by_frame: list[list[tuple[int, str, str]]] = []
    for frame in range(3):
        frame_codons: list[tuple[int, str, str]] = []
        for start_index in range(frame, len(sequence) - 2, 3):
            codon = sequence[start_index : start_index + 3]
            if codon == START_CODON:
                frame_codons.append((start_index + 1, codon, "Start"))
            elif codon in STOP_CODONS:
                frame_codons.append((start_index + 1, codon, "Stop"))
        codons_by_frame.append(frame_codons)
        start_count = sum(kind == "Start" for _, _, kind in frame_codons)
        stop_count = sum(kind == "Stop" for _, _, kind in frame_codons)
        print(
            f"  {FRAME_NAMES[frame]}: {start_count:,} start codons, "
            f"{stop_count:,} stop codons"
        )
    return codons_by_frame


def create_figure(
    fasta_header: str,
    sequence: str,
    cds_start: int | None,
    cds_end: int | None,
    cds_frame: int | None,
    codons_by_frame: Sequence[Sequence[tuple[int, str, str]]],
) -> go.Figure:
    """Create the nucleotide, CDS, and three reading-frame Plotly tracks."""
    figure = go.Figure()
    sequence_length = len(sequence)
    positions = list(range(1, sequence_length + 1))

    # Place one hoverable letter at each transcript nucleotide position.
    figure.add_trace(go.Scatter(
        x=positions,
        y=[4.1] * sequence_length,
        text=list(sequence),
        mode="text",
        textfont=dict(
            color=[NUCLEOTIDE_COLORS.get(base, "#777777") for base in sequence],
            size=9,
        ),
        customdata=positions,
        hovertemplate="Position: %{customdata} nt<br>Nucleotide: %{text}<extra></extra>",
        showlegend=False,
    ))

    # Draw the three frame bands and the annotated CDS rectangle.
    for frame, (y_bottom, y_top) in enumerate(FRAME_BOUNDS):
        figure.add_shape(
            type="rect",
            x0=1,
            x1=sequence_length,
            y0=y_bottom,
            y1=y_top,
            fillcolor=FRAME_COLORS[frame],
            line=dict(width=0),
            layer="below",
        )

    if cds_start is not None and cds_end is not None and cds_frame is not None:
        cds_color = FRAME_COLORS[cds_frame]
        figure.add_shape(
            type="rect",
            x0=cds_start,
            x1=cds_end,
            y0=3.05,
            y1=3.75,
            fillcolor=cds_color,
            line=dict(color=cds_color, width=1),
            layer="below",
        )
        figure.add_trace(go.Scatter(
            x=[(cds_start + cds_end) / 2],
            y=[3.4],
            text=[fasta_header],
            mode="text",
            textfont=dict(color="black", size=11),
            hoverinfo="skip",
            showlegend=False,
        ))

    # Add visible vertical codon markers and transparent hover targets.
    for frame, frame_codons in enumerate(codons_by_frame):
        y_bottom, y_top = FRAME_BOUNDS[frame]
        hover_positions: list[int] = []
        hover_text: list[str] = []
        for position, codon, codon_type in frame_codons:
            marker_color = "white" if codon_type == "Start" else "black"
            figure.add_shape(
                type="line",
                x0=position,
                x1=position,
                y0=y_bottom,
                y1=y_top,
                line=dict(color=marker_color, width=2),
                layer="above",
            )
            hover_positions.append(position)
            hover_text.append(
                f"Position: {position} nt<br>Codon: {codon}<br>"
                f"Type: {codon_type}<br>Frame: {FRAME_NAMES[frame]}"
            )

        if hover_positions:
            figure.add_trace(go.Scatter(
                x=hover_positions,
                y=[(y_bottom + y_top) / 2] * len(hover_positions),
                mode="markers",
                marker=dict(color="rgba(0,0,0,0.01)", size=12),
                text=hover_text,
                hovertemplate="%{text}<extra></extra>",
                showlegend=False,
            ))

    figure.update_layout(
        width=1200,
        height=200,
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(l=20, r=20, t=10, b=45),
        xaxis=dict(
            range=[1, sequence_length],
            title="position [nt]",
            showgrid=False,
            zeroline=False,
            ticks="outside",
            showline=False,
        ),
        yaxis=dict(
            range=[-0.1, 4.5],
            visible=False,
            fixedrange=True,
        ),
        showlegend=False,
        hovermode="closest",
    )
    return figure


def save_html(figure: go.Figure, fasta_header: str, output_path: Path) -> None:
    """Write standalone HTML with the FASTA header above the embedded plot."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing standalone HTML: {output_path}")
    plot_html = figure.to_html(include_plotlyjs=True, full_html=False)
    safe_header = html_escape(fasta_header)
    html_document = (
        "<!doctype html>\n"
        "<html><head><meta charset=\"utf-8\"></head><body>\n"
        f"<h2>{safe_header}</h2>\n"
        f"{plot_html}\n"
        "</body></html>\n"
    )
    output_path.write_text(html_document, encoding="utf-8")


def main() -> None:
    arguments = parse_args()
    fasta_header, sequence = read_fasta(arguments.fasta)
    transcript_id = fasta_header.split(maxsplit=1)[0].split("|")[-1]
    cds_start = arguments.cds_start
    cds_end = arguments.cds_end
    cds_length: int | None = None
    cds_frame: int | None = None

    if (cds_start is None) != (cds_end is None):
        print(
            "Error: provide both --cds_start and --cds_end, or omit both; "
            "continuing without an annotated CDS rectangle.",
            file=sys.stderr,
        )
        cds_start = None
        cds_end = None
    elif cds_start is None and cds_end is None:
        coordinates = find_cds_coordinates(transcript_id)
        if coordinates is None:
            print(
                f"Error: CDS coordinates for {transcript_id} were not found; "
                "continuing without an annotated CDS rectangle.",
                file=sys.stderr,
            )
        else:
            cds_start, cds_end = coordinates

    if cds_start is not None and cds_end is not None:
        try:
            cds_length, cds_frame = validate_input(sequence, cds_start, cds_end)
        except ValueError as error:
            print(
                f"Error: {error} Continuing without an annotated CDS rectangle.",
                file=sys.stderr,
            )
            cds_start = None
            cds_end = None
            cds_length = None
            cds_frame = None
        else:
            print(
                f"Validated CDS: {cds_start}–{cds_end} nt "
                f"({cds_length:,} nt), Frame +{cds_frame}"
            )

    codons_by_frame = identify_codons(sequence)
    figure = create_figure(
        fasta_header,
        sequence,
        cds_start,
        cds_end,
        cds_frame,
        codons_by_frame,
    )
    save_html(figure, fasta_header, arguments.output)

    counts = [
        Counter(codon_type for _, _, codon_type in frame_codons)
        for frame_codons in codons_by_frame
    ]
    print("Summary:")
    print(f"  Transcript ID: {transcript_id}")
    print(f"  Sequence length: {len(sequence):,} nt")
    if cds_start is not None and cds_end is not None and cds_length is not None and cds_frame is not None:
        print(f"  CDS coordinates: {cds_start}–{cds_end} nt")
        print(f"  CDS length: {cds_length:,} nt")
        print(f"  CDS reading frame: Frame +{cds_frame}")
    else:
        print("  CDS coordinates: unavailable; no annotated CDS rectangle")
    for frame, frame_counts in enumerate(counts):
        print(
            f"  {FRAME_NAMES[frame]}: {frame_counts['Start']:,} start, "
            f"{frame_counts['Stop']:,} stop codons"
        )
    print(f"  Output HTML: {arguments.output}")


if __name__ == "__main__":
    main()
