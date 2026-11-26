#!/usr/bin/env python3
"""Plot FASTA transcripts, annotated CDS regions, and codons in three frames."""

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
            "Create a standalone Plotly HTML page with vertically stacked plots "
            "for transcripts in a FASTA file."
        )
    )
    parser.add_argument("--fasta", required=True, type=Path, help="Input FASTA file with one or more transcripts")
    parser.add_argument("--output", required=True, type=Path, help="Output standalone HTML file")
    return parser.parse_args()


def read_fasta(fasta_path: Path) -> list[tuple[str, str]]:
    """Read all FASTA records and return each full header and sequence."""
    if not fasta_path.is_file():
        raise FileNotFoundError(f"FASTA file does not exist: {fasta_path}")

    header: str | None = None
    sequence_parts: list[str] = []
    records: list[tuple[str, str]] = []
    line_count = 0
    print(f"Reading FASTA: {fasta_path}")
    with fasta_path.open(encoding="utf-8") as fasta_file:
        for line_count, line in enumerate(fasta_file, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(sequence_parts).upper()
                    if not sequence:
                        raise ValueError(f"FASTA record {header!r} has an empty sequence.")
                    records.append((header, sequence))
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"FASTA header is empty at line {line_count}.")
                sequence_parts = []
            elif header is None:
                raise ValueError(f"Sequence data found before the FASTA header at line {line_count}.")
            else:
                sequence_parts.append("".join(line.split()))
            if line_count % 100_000 == 0:
                print(f"  Read {line_count:,} FASTA lines")

    if header is None:
        raise ValueError(f"No FASTA record found in {fasta_path}.")
    sequence = "".join(sequence_parts).upper()
    if not sequence:
        raise ValueError(f"FASTA record {header!r} has an empty sequence.")
    records.append((header, sequence))
    print(f"  Loaded {len(records):,} transcript records")
    return records


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


def find_cds_coordinates(
    transcript_ids: set[str],
) -> dict[str, tuple[int, int] | None]:
    """Search canonical transcript TSVs for the requested CDS intervals."""
    annotation_dir = Path("data/large_files/ensembl/gtf")
    annotation_files = sorted(annotation_dir.glob("*.canonical_transcripts.tsv"))
    if not annotation_files:
        print(
            f"Error: no *.canonical_transcripts.tsv files found in {annotation_dir}.",
            file=sys.stderr,
        )
        return {}

    print(
        f"Looking up CDS coordinates for {len(transcript_ids):,} transcripts "
        f"in {len(annotation_files):,} annotation TSV files"
    )
    coordinates: dict[str, tuple[int, int] | None] = {}
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
                transcript_id = row["transcript_id"].strip()
                if transcript_id not in transcript_ids or transcript_id in coordinates:
                    continue
                try:
                    start = int(row["cds_start"])
                    end = int(row["cds_end"])
                except (TypeError, ValueError):
                    coordinates[transcript_id] = None
                    continue
                if start <= 0 or end <= 0:
                    coordinates[transcript_id] = None
                    continue
                coordinates[transcript_id] = (start, end)
                print(f"  Found CDS coordinates {start}–{end} for {transcript_id} in {annotation_path.name}")

        if len(coordinates) == len(transcript_ids):
            break

    return coordinates


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


def save_html(plots: Sequence[tuple[str, go.Figure]], output_path: Path) -> None:
    """Write one standalone HTML page with a vertically stacked plot per record."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing standalone HTML: {output_path}")
    sections: list[str] = []
    for plot_index, (fasta_header, figure) in enumerate(plots):
        plot_html = figure.to_html(include_plotlyjs=(plot_index == 0), full_html=False)
        safe_header = html_escape(fasta_header)
        sections.append(
            f'<section style="margin: 0 0 32px 0;">\n'
            f"<h2>{safe_header}</h2>\n{plot_html}\n</section>"
        )
    html_document = "\n".join((
        "<!doctype html>",
        "<html><head><meta charset=\"utf-8\"></head><body>",
        *sections,
        "</body></html>",
        "",
    ))
    output_path.write_text(html_document, encoding="utf-8")


def main() -> None:
    arguments = parse_args()
    records = read_fasta(arguments.fasta)
    transcript_ids = {
        header.split(maxsplit=1)[0].split("|")[-1]
        for header, _ in records
    }
    cds_by_transcript = find_cds_coordinates(transcript_ids)
    plots: list[tuple[str, go.Figure]] = []

    for record_index, (fasta_header, sequence) in enumerate(records, start=1):
        transcript_id = fasta_header.split(maxsplit=1)[0].split("|")[-1]
        print(f"Processing transcript {record_index:,}/{len(records):,}: {transcript_id}")
        coordinates = cds_by_transcript.get(transcript_id)
        cds_start: int | None = None
        cds_end: int | None = None
        cds_length: int | None = None
        cds_frame: int | None = None

        if coordinates is None:
            print(
                f"Error: CDS coordinates for {transcript_id} were not found or are incomplete; "
                "plotting without an annotated CDS rectangle.",
                file=sys.stderr,
            )
        else:
            cds_start, cds_end = coordinates
            try:
                cds_length, cds_frame = validate_input(sequence, cds_start, cds_end)
            except ValueError as error:
                print(
                    f"Error: {error} Plotting without an annotated CDS rectangle.",
                    file=sys.stderr,
                )
                cds_start = None
                cds_end = None
                cds_length = None
                cds_frame = None
            else:
                print(
                    f"  CDS {cds_start}–{cds_end} nt ({cds_length:,} nt), "
                    f"Frame +{cds_frame}"
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
        plots.append((fasta_header, figure))

        for frame, frame_codons in enumerate(codons_by_frame):
            counts = Counter(codon_type for _, _, codon_type in frame_codons)
            print(
                f"  {FRAME_NAMES[frame]}: {counts['Start']:,} start, "
                f"{counts['Stop']:,} stop codons"
            )

    save_html(plots, arguments.output)
    print(f"Wrote {len(plots):,} transcript plots to {arguments.output}")


if __name__ == "__main__":
    main()
