#!/usr/bin/env python3
"""Plot GC content against predicted RNA structure strength for human transcripts."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import pandas as pd

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as error:
    raise ImportError(
        "This script requires matplotlib, which is not currently listed in scripts-env.yml."
    ) from error


ROOT = Path(__file__).resolve().parents[2]
INPUT_FILE = ROOT / "results" / "human" / "05_mane_start_codon_structure.tsv"
OUTPUT_FILE = ROOT / "results" / "human" / "06_gc_vs_mfe.png"
REQUIRED_COLUMNS = ("gene_name", "transcript", "GC", "free_energy")
GABARAP_TRANSCRIPT = "NM_007278.2"


def load_plot_data(path: Path) -> tuple[pd.DataFrame, int]:
    """Load results, coerce numeric columns, and count rows excluded as invalid."""
    if not path.is_file():
        raise FileNotFoundError(f"Step 5 results not found: {path}")
    frame = pd.read_csv(path, sep="\t", dtype={"gene_name": str, "transcript": str})
    missing_columns = set(REQUIRED_COLUMNS) - set(frame.columns)
    if missing_columns:
        raise ValueError(f"Input is missing required columns: {sorted(missing_columns)}")

    frame["GC"] = pd.to_numeric(frame["GC"], errors="coerce")
    frame["free_energy"] = pd.to_numeric(frame["free_energy"], errors="coerce")
    valid = (
        frame["GC"].notna()
        & frame["free_energy"].notna()
        & frame["transcript"].fillna("").str.strip().ne("")
    )
    invalid_count = int((~valid).sum())
    frame = frame.loc[valid].copy()
    frame["structure_strength"] = -frame["free_energy"]
    frame["is_gabarap"] = (
        frame["gene_name"].fillna("").str.strip().str.casefold().eq("gabarap")
        | frame["transcript"].eq(GABARAP_TRANSCRIPT)
    )
    return frame, invalid_count


def plot_data(frame: pd.DataFrame, output_path: Path) -> None:
    """Create and atomically save a publication-style 300-DPI scatter plot."""
    if frame.empty:
        raise ValueError("No rows with valid transcript, GC, and free_energy values to plot")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    other = frame.loc[~frame["is_gabarap"]]
    gabarap = frame.loc[frame["is_gabarap"]]

    fig, ax = plt.subplots(figsize=(10, 7), layout="constrained")
    ax.scatter(
        other["GC"],
        other["structure_strength"],
        s=24,
        color="#718096",
        alpha=0.58,
        edgecolors="none",
        label="Other MANE Select transcripts",
        rasterized=True,
    )
    if not gabarap.empty:
        ax.scatter(
            gabarap["GC"],
            gabarap["structure_strength"],
            s=75,
            color="#C43C39",
            edgecolors="white",
            linewidths=0.8,
            label="GABARAP",
            zorder=3,
        )
        for _, row in gabarap.iterrows():
            ax.annotate(
                "GABARAP",
                (row["GC"], row["structure_strength"]),
                xytext=(8, 8),
                textcoords="offset points",
                color="#A52A2A",
                fontsize=10,
                fontweight="bold",
            )

    ax.set_xlabel("GC content (%)")
    ax.set_ylabel("Predicted RNA structure strength, −MFE (kcal/mol)")
    ax.set_title("GC Content versus Predicted RNA Structure Strength\nHuman MANE Select Transcripts")
    ax.grid(True, color="#D9DEE5", linewidth=0.7, alpha=0.8)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, loc="best")

    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            suffix=".png", dir=output_path.parent, delete=False
        ) as temporary:
            temporary_path = Path(temporary.name)
        fig.savefig(temporary_path, dpi=300, format="png", facecolor="white")
        os.replace(temporary_path, output_path)
    finally:
        plt.close(fig)
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def main() -> int:
    try:
        frame, invalid_count = load_plot_data(INPUT_FILE)
        plot_data(frame, OUTPUT_FILE)
    except (OSError, ValueError, pd.errors.ParserError) as error:
        print(f"Plot creation failed: {error}")
        return 1

    print(f"Transcripts plotted: {len(frame)}")
    print(f"Rows skipped for missing or non-numeric values: {invalid_count}")
    gabarap = frame.loc[frame["is_gabarap"]]
    if gabarap.empty:
        print("GABARAP transcript was not found among valid rows.")
    else:
        for _, row in gabarap.iterrows():
            print(
                f"GABARAP ({row['transcript']}): GC={row['GC']:.1f}%, "
                f"-MFE={row['structure_strength']:.2f} kcal/mol"
            )
    print(f"Plot saved to: {OUTPUT_FILE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
