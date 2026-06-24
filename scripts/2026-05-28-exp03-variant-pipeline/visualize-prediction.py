#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 12:00:29 2026

@author: e6260
"""
"""
Visualize AlphaGenome predict_variant outputs.

This script loads a previously saved AlphaGenome prediction object
and generates local biological visualizations for:
- RNA-seq
- splice sites
- splice site usage
- splice junctions

The script is intended to complement run-variant-pipeline.py.
"""

import argparse
import pickle
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

REF_COLOR = "#1f77b4"      # azul
ALT_COLOR = "#7a0019"      # granate
DELTA_COLOR = "#ff8c00"    # naranja
EXON_COLOR = "#1f77b4"     # azul
INTRON_COLOR = "#7a0019"   # granate
VARIANT_COLOR = "#D4A017"  # amarillo intenso (gold)

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Visualize AlphaGenome predict_variant outputs from prediction.pkl"
    )

    parser.add_argument(
        "--prediction-pkl",
        required=True,
        help="Path to saved AlphaGenome prediction.pkl file"
    )

    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory where figures will be saved. Default: <prediction_dir>/figures"
    )

    parser.add_argument(
        "--runlog",
        default=None,
        help="Optional path to runlog.txt. Default: same directory as prediction.pkl"
    )

    parser.add_argument(
        "--region-start",
        type=int,
        default=None,
        help="Optional local region start, 1-based"
    )

    parser.add_argument(
        "--region-end",
        type=int,
        default=None,
        help="Optional local region end, 1-based"
    )

    parser.add_argument(
        "--view-flank",
        type=int,
        default=None,
        help=(
            "Optional visualization flank around the variant position, in bp. "
            "If provided, it overrides the runlog local window for plotting only."
        )
    )
    parser.add_argument(
        "--variant-pos",
        type=int,
        default=None,
        help="Optional variant position, 1-based. If omitted, read from runlog."
    )

    parser.add_argument(
        "--strand",
        default="both",
        choices=["+", "-", "both"],
        help="Strand to visualize: +, -, or both"
    )

    parser.add_argument(
        "--title",
        default="AlphaGenome local prediction",
        help="Figure title prefix"
    )

    return parser.parse_args()


def load_prediction(prediction_pkl):
    with open(prediction_pkl, "rb") as f:
        prediction = pickle.load(f)

    return prediction


def parse_runlog_value(runlog_path, key):
    """
    Extract a value from runlog lines formatted as:
    key: value
    """
    if runlog_path is None or not runlog_path.exists():
        return None

    with open(runlog_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(f"{key}:"):
                return line.split(":", 1)[1].strip()

    return None


def resolve_region(args, prediction_dir):
    """
    Resolve visualization region.

    Priority:
    1. Manual --region-start and --region-end
    2. local_start_1based and local_end_1based from runlog.txt
    3. fallback around variant_pos
    """

    runlog_path = Path(args.runlog) if args.runlog else prediction_dir / "runlog.txt"

    if args.region_start is not None and args.region_end is not None:
        region_start = args.region_start
        region_end = args.region_end
        region_source = "manual_arguments"

    else:
        local_start = parse_runlog_value(runlog_path, "local_start_1based")
        local_end = parse_runlog_value(runlog_path, "local_end_1based")

        if local_start is not None and local_end is not None:
            region_start = int(local_start)
            region_end = int(local_end)
            region_source = "runlog_local_window"

        else:
            variant_pos = args.variant_pos
            if variant_pos is None:
                variant_pos_runlog = parse_runlog_value(runlog_path, "position_1based")
                if variant_pos_runlog is not None:
                    variant_pos = int(variant_pos_runlog)

            if variant_pos is None:
                raise ValueError(
                    "Could not resolve region. Provide --region-start and --region-end, "
                    "or provide a valid runlog.txt, or provide --variant-pos."
                )

            flank = 5000
            region_start = max(1, variant_pos - flank)
            region_end = variant_pos + flank
            region_source = "fallback_variant_pos_plus_minus_5kb"
            
    if args.view_flank is not None:
        variant_pos = args.variant_pos
        
        if variant_pos is None:
            variant_pos_runlog = parse_runlog_value(runlog_path, "position_1based")
            
            if variant_pos_runlog is not None:
                variant_pos = int(variant_pos_runlog)

        if variant_pos is None:
            raise ValueError(
                "--view-flank was provided, but variant position could not be resolved."
        )

        region_start = max(1, variant_pos - args.view_flank)
        region_end = variant_pos + args.view_flank
        region_source = f"variant_pos_plus_minus_{args.view_flank}bp"

    return region_start, region_end, region_source, runlog_path


def resolve_variant_position(args, runlog_path):
    if args.variant_pos is not None:
        return args.variant_pos

    value = parse_runlog_value(runlog_path, "position_1based")
    if value is not None:
        return int(value)

    return None

def trackdata_to_local_dataframe(ref_data, alt_data, region_start, region_end, track_index=0):
    """
    Convert one AlphaGenome TrackData track into a local REF/ALT dataframe.

    AlphaGenome stores predictions as arrays indexed relative to the prediction interval.
    This function converts array indices into genomic 1-based coordinates and filters
    the selected local region.
    """

    if ref_data is None or alt_data is None:
        raise ValueError("REF or ALT TrackData is missing.")

    if track_index >= ref_data.values.shape[1]:
        raise ValueError(
            f"track_index {track_index} is not available. "
            f"Available tracks: 0-{ref_data.values.shape[1] - 1}"
        )

    n_positions = ref_data.values.shape[0]

    genomic_position_1based = (
        ref_data.interval.start
        + pd.Series(range(n_positions))
        + 1
    )

    df = pd.DataFrame({
        "genomic_position_1based": genomic_position_1based,
        "ref_value": ref_data.values[:, track_index],
        "alt_value": alt_data.values[:, track_index],
    })

    df["delta_alt_ref"] = df["alt_value"] - df["ref_value"]

    local_df = df[
        (df["genomic_position_1based"] >= region_start)
        &
        (df["genomic_position_1based"] <= region_end)
    ].copy()

    return local_df


def plot_rna_seq_ref_alt_and_delta(
    prediction,
    region_start,
    region_end,
    variant_pos,
    output_dir,
    title,
    annotation_df=None
):
    """
    Plot RNA-seq REF vs ALT, ALT-REF delta, and optional gene annotation track.
    """

    ref_data = prediction.reference.rna_seq
    alt_data = prediction.alternate.rna_seq

    if ref_data is None or alt_data is None:
        print("Skipping RNA-seq: output not available.")
        return None

    df = trackdata_to_local_dataframe(
        ref_data=ref_data,
        alt_data=alt_data,
        region_start=region_start,
        region_end=region_end,
        track_index=0
    )

    if df.empty:
        print("Skipping RNA-seq: no data in selected local region.")
        return None

    track_name = ref_data.metadata.iloc[0].get("name", "RNA-seq")

    if annotation_df is not None:
        fig, axes = plt.subplots(
            nrows=3,
            ncols=1,
            figsize=(12, 7),
            sharex=True,
            gridspec_kw={"height_ratios": [2, 1, 0.6]}
        )
        ax_signal, ax_delta, ax_annot = axes
    else:
        fig, axes = plt.subplots(
            nrows=2,
            ncols=1,
            figsize=(12, 6),
            sharex=True,
            gridspec_kw={"height_ratios": [2, 1]}
        )
        ax_signal, ax_delta = axes
        ax_annot = None

    ax_signal.plot(df["genomic_position_1based"], df["ref_value"], label="REF", color=REF_COLOR)
    ax_signal.plot(df["genomic_position_1based"], df["alt_value"], label="ALT", color=ALT_COLOR)
    ax_signal.set_ylabel("Predicted RNA-seq signal")
    ax_signal.set_title(f"{title} — RNA-seq REF vs ALT\n{track_name}")
    ax_signal.legend()

    ax_delta.axhline(0, linestyle="--", linewidth=1)
    ax_delta.fill_between(df["genomic_position_1based"], 0, df["delta_alt_ref"], color=DELTA_COLOR,
    alpha=0.35,)

    ax_delta.plot(df["genomic_position_1based"], df["delta_alt_ref"], color=DELTA_COLOR, linewidth=1, label="ALT - REF",)

    ax_delta.legend(loc="upper right")

    if variant_pos is not None:
        ax_signal.axvline(variant_pos, linestyle="-.", linewidth=1, color=VARIANT_COLOR, )
        ax_signal.text(variant_pos + 150, ax_signal.get_ylim()[1] * 0.95, f"KI\n({variant_pos:,})", color=VARIANT_COLOR, ha="left", va="top", fontweight="bold", fontsize=9, )
        ax_delta.axvline(variant_pos, linestyle="-.", linewidth=1, color=VARIANT_COLOR, )

    if ax_annot is not None:
        local_annotation = filter_annotation_to_region(
            annotation_df,
            region_start,
            region_end
        )

        ax_annot.set_ylim(-0.5, 1.8)
        ax_annot.set_yticks([0, 1])
        ax_annot.set_yticklabels(["introns", "exons"])
        ax_annot.tick_params(axis="x", bottom=False, labelbottom=False)
        ax_annot.set_xlabel("")
        ax_annot.set_title("SEC23B MANE Select annotation")

        for _, row in local_annotation.iterrows():
            start = row["plot_start"]
            width = row["plot_end"] - row["plot_start"] + 1

            if row["type"] == "exon":
                ax_annot.broken_barh(
                    [(start, width)],
                    (0.75, 0.5),
                    facecolors=EXON_COLOR,
                    edgecolors=EXON_COLOR,
                )   
              
                ax_annot.text(
                    start + width / 2,
                    1.35,
                    row["name"].replace("exon ", "E"),
                    ha="center",
                    va="bottom",
                    fontsize=8,
                    color="black",
                )   
            elif row["type"] == "intron":
                ax_annot.hlines(
                    y=0,
                    xmin=start,
                    xmax=row["plot_end"],
                    linewidth=1.5,
                    color=INTRON_COLOR,
                )
                ax_annot.text(
                    start + width / 2,
                    0.15,
                    row["name"].replace("intron ", "I"),
                    ha="center",
                    va="bottom",
                    fontsize=7,
                    color=INTRON_COLOR
                )

        if variant_pos is not None:
            ax_annot.axvline(variant_pos, linestyle="-.", linewidth=1, color=VARIANT_COLOR,)

    else:
        ax_delta.set_xlabel("Genomic position (1-based)")

    fig.tight_layout()

    output_path = output_dir / "01_rna_seq_ref_alt_delta_annotation.png"
    fig.savefig(output_path, dpi=300)
    plt.close(fig)

    print(f"RNA-seq figure saved to: {output_path}")

    return output_path

def get_sec23b_mane_annotation():
    """
    Return SEC23B MANE Select exon/intron coordinates.
    Coordinates are GRCh38, 1-based.
    """

    regions = [
        {"name": "exon 1", "start": 18507940, "end": 18507972, "type": "exon"},
        {"name": "intron 1-2", "start": 18507973, "end": 18510821, "type": "intron"},
        {"name": "exon 2", "start": 18510822, "end": 18511056, "type": "exon"},
        {"name": "intron 2-3", "start": 18511057, "end": 18512224, "type": "intron"},
        {"name": "exon 3", "start": 18512225, "end": 18512282, "type": "exon"},
        {"name": "intron 3-4", "start": 18512283, "end": 18515649, "type": "intron"},
        {"name": "exon 4", "start": 18515650, "end": 18515736, "type": "exon"},
        {"name": "intron 4-5", "start": 18515737, "end": 18524432, "type": "intron"},
        {"name": "exon 5", "start": 18524433, "end": 18524669, "type": "exon"},
    ]

    return pd.DataFrame(regions)


def filter_annotation_to_region(annotation_df, region_start, region_end):
    """
    Keep only exon/intron annotations overlapping the plotted local region.
    """

    local_annotation = annotation_df[
        (annotation_df["start"] <= region_end)
        &
        (annotation_df["end"] >= region_start)
    ].copy()

    local_annotation["plot_start"] = local_annotation["start"].clip(lower=region_start)
    local_annotation["plot_end"] = local_annotation["end"].clip(upper=region_end)

    return local_annotation

def main():
    args = parse_arguments()

    prediction_pkl = Path(args.prediction_pkl)
    prediction_dir = prediction_pkl.parent

    output_dir = Path(args.output_dir) if args.output_dir else prediction_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)

    prediction = load_prediction(prediction_pkl)

    region_start, region_end, region_source, runlog_path = resolve_region(
        args,
        prediction_dir
    )

    variant_pos = resolve_variant_position(args, runlog_path)

    print("Prediction loaded successfully.")
    print(f"Prediction file: {prediction_pkl}")
    print(f"Output directory: {output_dir}")
    print(f"Runlog: {runlog_path}")
    print(f"Region: {region_start}-{region_end}")
    print(f"Region source: {region_source}")
    print(f"Variant position: {variant_pos}")
    print(f"Strand: {args.strand}")

    print("\nAvailable reference outputs:")
    for output_name in [
        "rna_seq",
        "splice_sites",
        "splice_site_usage",
        "splice_junctions",
    ]:
        output = getattr(prediction.reference, output_name, None)
        if output is None:
            print(f"- {output_name}: not available")
        else:
            print(f"- {output_name}: available")
            
    sec23b_annotation = get_sec23b_mane_annotation()

    plot_rna_seq_ref_alt_and_delta(
        prediction=prediction,
        region_start=region_start,
        region_end=region_end,
        variant_pos=variant_pos,
        output_dir=output_dir,
        title=args.title,
        annotation_df=sec23b_annotation
    )
if __name__ == "__main__":
    main()