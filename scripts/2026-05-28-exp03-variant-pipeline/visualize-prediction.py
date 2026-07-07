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
from datetime import datetime

import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from matplotlib.lines import Line2D


REF_COLOR = "#1f77b4"      # azul
ALT_COLOR = "#7a0019"      # granate
DELTA_COLOR = "#ff8c00"    # naranja
DELTA_SECOND_COLOR = REF_COLOR

EXON_COLOR = "#1f77b4"     # azul
INTRON_COLOR = "#7a0019"   # granate
VARIANT_COLOR = "#D4A017"  # amarillo intenso (gold)

RNA_SEQ_LABEL = "RNA-seq"
SPLICE_SITES_LABEL = "Splice sites"
SPLICE_SITE_USAGE_LABEL = "Splice site usage"
SPLICE_JUNCTIONS_LABEL = "Splice junctions"


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
    "--annotation-gtf",
    default=None,
    help=(
        "Optional Ensembl/GENCODE GTF file used to draw exon/intron "
        "annotation below the plots. Must match the genome assembly used "
        "by AlphaGenome, normally GRCh38."
        )
    )

    parser.add_argument(
    "--gene-name",
    default=None,
    help=(
        "Optional gene name used to extract exon/intron annotation "
        "from the GTF, e.g. SEC23B or DLG1."
        )
    )

    parser.add_argument(
    "--transcript-id",
    default=None,
    help=(
        "Optional transcript ID used to select a specific transcript "
        "from the GTF. If omitted, the script will use the best available "
        "transcript for the selected gene."
        )
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
        "--junction-min-positive-delta",
        type=float,
        default=0.5,
        help="Minimum positive ALT-REF delta required to display increased junctions."
    )

    parser.add_argument(
        "--junction-max-negative-delta",
        type=float,
        default=-0.5,
        help="Maximum negative ALT-REF delta required to display decreased junctions."
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


def validate_existing_file(path, argument_name):
    """
    Validate that an input file exists.
    """

    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(
            f"{argument_name} file was not found: {path}"
        )

    if not path.is_file():
        raise ValueError(
            f"{argument_name} must point to a file, not a directory: {path}"
        )

    return path

def resolve_runlog_path(args, prediction_dir):
    """
    Resolve and validate the runlog path.

    The runlog is optional only if region and variant position are provided
    manually by the user.
    """

    if args.runlog is not None:
        return validate_existing_file(
            args.runlog,
            "--runlog"
        )

    default_runlog_path = prediction_dir / "runlog.txt"

    if default_runlog_path.exists():
        return default_runlog_path

    return None

def validate_region(region_start, region_end):
    """
    Validate plotting region coordinates.

    Coordinates must be positive 1-based integers and start must be <= end.
    """

    if region_start < 1 or region_end < 1:
        raise ValueError(
            "Plotting region coordinates must be positive 1-based integers."
        )

    if region_start > region_end:
        raise ValueError(
            f"Invalid plotting region: start ({region_start}) is greater "
            f"than end ({region_end})."
        )

    return region_start, region_end


def validate_variant_inside_region(variant_pos, region_start, region_end):
    """
    Validate that the variant position is inside the plotted region.
    """

    if variant_pos is None:
        return

    if not (region_start <= variant_pos <= region_end):
        raise ValueError(
            f"Variant position ({variant_pos}) is outside the plotted region "
            f"({region_start}-{region_end}). Adjust --region-start/--region-end "
            "or use --view-flank."
        )

def validate_annotation_arguments(args):
    """
    Validate gene annotation arguments.

    Gene annotation from GTF requires both:
    - --annotation-gtf
    - --gene-name

    If neither is provided, figures are generated without exon/intron annotation.
    """

    if args.annotation_gtf is None and args.gene_name is None:
        return

    if args.annotation_gtf is None and args.gene_name is not None:
        raise ValueError(
            "--gene-name was provided, but --annotation-gtf is missing. "
            "Provide a GTF file, for example: "
            "--annotation-gtf raw-data/Homo_sapiens.GRCh38.115.gtf"
        )

    if args.annotation_gtf is not None and args.gene_name is None:
        raise ValueError(
            "--annotation-gtf was provided, but --gene-name is missing. "
            "Provide the gene name to extract, for example: --gene-name SEC23B."
        )

    args.annotation_gtf = validate_existing_file(
        args.annotation_gtf,
        "--annotation-gtf"
    )

def get_available_prediction_outputs(prediction):
    """
    Check availability of the AlphaGenome predict_variant outputs used by
    this visualization script.
    """

    output_names = [
        "rna_seq",
        "splice_sites",
        "splice_site_usage",
        "splice_junctions",
    ]

    available_outputs = {}

    for output_name in output_names:
        ref_data = getattr(prediction.reference, output_name, None)
        alt_data = getattr(prediction.alternate, output_name, None)

        available_outputs[output_name] = (
            ref_data is not None
            and alt_data is not None
        )

    return available_outputs


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


def resolve_region(args, runlog_path):
    """
    Resolve visualization region.

    Priority:
    1. Manual --region-start and --region-end
    2. local_start_1based and local_end_1based from runlog.txt
    3. fallback around variant_pos
    """

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
    region_start, region_end = validate_region(
        region_start,
        region_end
    )

    return region_start, region_end, region_source


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

def plot_mane_annotation_track(
    ax,
    annotation_df,
    region_start,
    region_end,
    variant_pos=None,
    transcript_info=None,
):
    """
    Plot MANE Select exon/intron annotation in a dedicated axis.
    """

    local_annotation = filter_annotation_to_region(
        annotation_df,
        region_start,
        region_end
    )

    ax.set_ylim(-0.5, 1.8)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["introns", "exons"])
    if transcript_info is not None:
        gene_name = transcript_info.get("gene_name", "Gene")
        selection_method = transcript_info.get("selection_method", "selected_transcript")
        selection_label = selection_method.replace("_", " ")

        ax.set_title(f"{gene_name} {selection_label} annotation")
    else:
        ax.set_title("Gene annotation")

    if transcript_info is not None:
        gene_strand = transcript_info.get("strand")

        if gene_strand == "-":
            ax.text(
                0.99,
                1.05,
                "Gene strand: − | transcription direction: ←",
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=8,
                color="black",
            )

        elif gene_strand == "+":
            ax.text(
                0.99,
                1.05,
                "Gene strand: + | transcription direction: →",
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=8,
                color="black",
            )


    ax.tick_params(
        axis="x",
        bottom=False,
        labelbottom=False
    )

    ax.set_xlabel("")

    for _, row in local_annotation.iterrows():
        start = row["plot_start"]
        width = row["plot_end"] - row["plot_start"] + 1

        if row["type"] == "exon":
            ax.broken_barh(
                [(start, width)],
                (0.75, 0.5),
                facecolors=EXON_COLOR,
                edgecolors=EXON_COLOR,
            )

            ax.text(
                start + width / 2,
                1.35,
                row["name"].replace("exon ", "E"),
                ha="center",
                va="bottom",
                fontsize=8,
                color="black",
            )

        elif row["type"] == "intron":
            ax.hlines(
                y=0,
                xmin=start,
                xmax=row["plot_end"],
                linewidth=1.5,
                color=INTRON_COLOR,
            )

            ax.text(
                start + width / 2,
                0.15,
                row["name"].replace("intron ", "I"),
                ha="center",
                va="bottom",
                fontsize=7,
                color=INTRON_COLOR,
            )

    if variant_pos is not None:
        ax.axvline(
            variant_pos,
            linestyle="-.",
            linewidth=1,
            color=VARIANT_COLOR,
        )

def get_splice_site_track_indices(strand):
    """
    Return donor and acceptor track indices for AlphaGenome splice_sites.

    AlphaGenome splice_sites metadata:
    track 0: donor +
    track 1: acceptor +
    track 2: donor -
    track 3: acceptor -
    """

    if strand == "+":
        return {
            "donor": 0,
            "acceptor": 1,
            "strand_label": "plus strand"
        }

    if strand == "-":
        return {
            "donor": 2,
            "acceptor": 3,
            "strand_label": "minus strand"
        }

    raise ValueError(
        "splice_sites plotting requires strand '+' or '-'. "
        "Use '--strand +' or '--strand -'."
    )
    
def plot_rna_seq_ref_alt_and_delta(
    prediction,
    region_start,
    region_end,
    variant_pos,
    output_dir,
    title,
    annotation_df=None,
    transcript_info=None
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
            gridspec_kw={"height_ratios": [0.6, 2, 1]}
        )
        ax_annot, ax_signal, ax_delta = axes
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
    ax_signal.set_title(f"{title} — {RNA_SEQ_LABEL}: REF vs ALT\n{track_name}")
    ax_signal.legend()

    ax_delta.axhline(0, linestyle="--", linewidth=1)
    ax_delta.fill_between(df["genomic_position_1based"], 0, df["delta_alt_ref"], color=DELTA_COLOR,
    alpha=0.35,)

    ax_delta.plot(df["genomic_position_1based"], df["delta_alt_ref"], color=DELTA_COLOR, linewidth=1, label="ALT - REF",)

    ax_delta.legend(loc="upper right")

    if variant_pos is not None:
        ax_signal.axvline(variant_pos, linestyle="-.", linewidth=1, color=VARIANT_COLOR, )
        ax_signal.text(variant_pos -150, ax_signal.get_ylim()[1] * 0.95, f"KI\n({variant_pos:,})", color=VARIANT_COLOR, ha="right", va="top", fontweight="bold", fontsize=9, )
        ax_delta.axvline(variant_pos, linestyle="-.", linewidth=1, color=VARIANT_COLOR, )

    if ax_annot is not None:
        plot_mane_annotation_track(
            ax=ax_annot,
            annotation_df=annotation_df,
            region_start=region_start,
            region_end=region_end,
            variant_pos=variant_pos,
            transcript_info=transcript_info,
        )

    else:
        ax_delta.set_xlabel("Genomic position (1-based)")

    fig.tight_layout()

    output_path = output_dir / "01_rna_seq.png"
    fig.savefig(output_path, dpi=300)
    plt.close(fig)

    print(f"{RNA_SEQ_LABEL} figure saved to: {output_path}")

    return output_path

def plot_splice_sites(
    prediction,
    region_start,
    region_end,
    variant_pos,
    output_dir,
    title,
    strand,
    annotation_df=None,
    transcript_info=None
):
    """
    Plot AlphaGenome splice_sites for donor and acceptor tracks
    in the selected strand.
    """

    ref_data = prediction.reference.splice_sites
    alt_data = prediction.alternate.splice_sites

    if ref_data is None or alt_data is None:
        print("Skipping splice_sites: output not available.")
        return None

    indices = get_splice_site_track_indices(strand)

    donor_df = trackdata_to_local_dataframe(
        ref_data=ref_data,
        alt_data=alt_data,
        region_start=region_start,
        region_end=region_end,
        track_index=indices["donor"]
    )

    acceptor_df = trackdata_to_local_dataframe(
        ref_data=ref_data,
        alt_data=alt_data,
        region_start=region_start,
        region_end=region_end,
        track_index=indices["acceptor"]
    )

    if donor_df.empty and acceptor_df.empty:
        print("Skipping splice_sites: no data in selected local region.")
        return None

    if annotation_df is not None:
        fig, axes = plt.subplots(
            nrows=4,
            ncols=1,
            figsize=(12, 8),
            sharex=True,
            gridspec_kw={"height_ratios": [0.6, 1.5, 1.5, 1]}
        )
        ax_annot, ax_donor, ax_acceptor, ax_delta = axes
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

    x_donor = donor_df["genomic_position_1based"]
    x_acceptor = acceptor_df["genomic_position_1based"]

    ax_donor.plot(x_donor, donor_df["ref_value"], label="REF donor", color=REF_COLOR, linestyle="-", linewidth=1.2,)
    ax_donor.plot(x_donor, donor_df["alt_value"], label="ALT donor", color=ALT_COLOR, linestyle="--", linewidth=1.2,)
    ax_donor.set_ylabel("Predicted donor probability")
    ax_donor.set_title(f"{title} — splice_sites REF vs ALT\n{indices['strand_label']}")
    ax_donor.legend(loc="lower right", bbox_to_anchor=(0.90,0.70),)
    
    ax_acceptor.plot(x_acceptor, acceptor_df["ref_value"], label="REF acceptor", color=REF_COLOR, linestyle="-", linewidth=1.2,)
    ax_acceptor.plot(x_acceptor, acceptor_df["alt_value"], label="ALT acceptor", color=ALT_COLOR, linestyle="--", linewidth=1.2,)
    ax_acceptor.set_ylabel("Predicted acceptor probability")
    ax_acceptor.legend(loc="lower right", bbox_to_anchor=(0.90, 0.70),)


    ax_delta.axhline(0, linestyle="--", linewidth=1)

    ax_delta.plot(
        x_donor,
        donor_df["delta_alt_ref"],
        color=DELTA_COLOR,
        linewidth=1,
        label="ALT - REF donor",
    )

    ax_delta.plot(
        x_acceptor,
        acceptor_df["delta_alt_ref"],
        color=DELTA_SECOND_COLOR,
        linewidth=1,
        linestyle="--",
        label="ALT - REF acceptor",
    )

    ax_delta.set_ylabel("ALT - REF")
    ax_delta.legend(loc="upper right")

    if variant_pos is not None:
        for ax in [ax_donor, ax_acceptor, ax_delta]:
            ax.axvline(
                variant_pos,
                linestyle="-.",
                linewidth=1,
                color=VARIANT_COLOR,
            )   

        ax_donor.text(
            variant_pos - 150,
            ax_donor.get_ylim()[1] * 0.95,
            f"KI\n({variant_pos:,})",
            color=VARIANT_COLOR,
            ha="right",
            va="top",
            fontweight="bold",
            fontsize=9,
        )     

    if ax_annot is not None:
        plot_mane_annotation_track(
            ax=ax_annot,
            annotation_df=annotation_df,
            region_start=region_start,
            region_end=region_end,
            variant_pos=variant_pos,
            transcript_info=transcript_info,
        )
    else:
        ax_delta.set_xlabel("Genomic position (1-based)") 

    fig.tight_layout()

    strand_suffix = "plus" if strand == "+" else "minus"
    output_path = output_dir / f"02_splice_sites_{strand_suffix}.png"

    fig.savefig(output_path, dpi=300)
    plt.close(fig)

    print(f"{SPLICE_SITES_LABEL} figure saved to: {output_path}")

    return output_path

def get_splice_site_usage_track_index(strand):
    """
    Return track index for AlphaGenome splice_site_usage.

    AlphaGenome splice_site_usage metadata:
    track 0: usage, plus strand
    track 1: usage, minus strand
    """

    if strand == "+":
        return {
            "track_index": 0,
            "strand_label": "plus strand"
        }

    if strand == "-":
        return {
            "track_index": 1,
            "strand_label": "minus strand"
        }

    raise ValueError(
        "splice_site_usage plotting requires strand '+' or '-'. "
        "Use '--strand +' or '--strand -'."
    )

def plot_splice_site_usage(
    prediction,
    region_start,
    region_end,
    variant_pos,
    output_dir,
    title,
    strand,
    annotation_df=None,
    transcript_info=None
):
    """
    Plot AlphaGenome splice_site_usage REF vs ALT and ALT-REF delta
    for the selected strand.
    """

    ref_data = prediction.reference.splice_site_usage
    alt_data = prediction.alternate.splice_site_usage

    if ref_data is None or alt_data is None:
        print("Skipping splice_site_usage: output not available.")
        return None

    usage_info = get_splice_site_usage_track_index(strand)

    df = trackdata_to_local_dataframe(
        ref_data=ref_data,
        alt_data=alt_data,
        region_start=region_start,
        region_end=region_end,
        track_index=usage_info["track_index"]
    )

    if df.empty:
        print("Skipping splice_site_usage: no data in selected local region.")
        return None

    if annotation_df is not None:
        fig, axes = plt.subplots(
            nrows=3,
            ncols=1,
            figsize=(12, 7),
            sharex=True,
            gridspec_kw={"height_ratios": [0.6, 2, 1]}
            )
        ax_annot, ax_signal, ax_delta = axes
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

    ax_signal.plot(
        df["genomic_position_1based"],
        df["ref_value"],
        label="REF",
        color=REF_COLOR,
        linewidth=1.2,
    )

    ax_signal.plot(
        df["genomic_position_1based"],
        df["alt_value"],
        label="ALT",
        color=ALT_COLOR,
        linestyle="--",
        linewidth=1.2,
    )

    ax_signal.set_ylabel("Predicted site usage")
    ax_signal.set_title(f"{title} — {SPLICE_SITE_USAGE_LABEL}: REF vs ALT\n{usage_info['strand_label']}")
    ax_signal.legend(loc="upper right", bbox_to_anchor=(0.90, 1.00))

    ax_delta.axhline(0, linestyle="--", linewidth=1)

    ax_delta.fill_between(
        df["genomic_position_1based"],
        0,
        df["delta_alt_ref"],
        color=DELTA_COLOR,
        alpha=0.35,
    )

    ax_delta.plot(
        df["genomic_position_1based"],
        df["delta_alt_ref"],
        color=DELTA_COLOR,
        linewidth=1,
        label="ALT - REF",
    )

    ax_delta.set_ylabel("ALT - REF")
    ax_delta.legend(loc="upper right")

    if variant_pos is not None:
        ax_signal.axvline(
            variant_pos,
            linestyle="-.",
            linewidth=1,
            color=VARIANT_COLOR,
        )

        ax_signal.text(
            variant_pos -150,
            ax_signal.get_ylim()[1] * 0.95,
            f"KI\n({variant_pos:,})",
            color=VARIANT_COLOR,
            ha="right",
            va="top",
            fontweight="bold",
            fontsize=9,
        )

        ax_delta.axvline(
            variant_pos,
            linestyle="-.",
            linewidth=1,
            color=VARIANT_COLOR,
        )

    if ax_annot is not None:
        plot_mane_annotation_track(
            ax=ax_annot,
            annotation_df=annotation_df,
            region_start=region_start,
            region_end=region_end,
            variant_pos=variant_pos,
            transcript_info=transcript_info,
        )
    else:
        ax_delta.set_xlabel("Genomic position (1-based)")

    fig.tight_layout()

    strand_suffix = "plus" if strand == "+" else "minus"
    output_path = output_dir / f"03_splice_site_usage_{strand_suffix}.png"

    fig.savefig(output_path, dpi=300)
    plt.close(fig)

    print(f"{SPLICE_SITE_USAGE_LABEL} figure saved to: {output_path}")

    return output_path

def parse_gtf_attributes(attribute_text):
    """
    Parse the attributes column of a GTF file into a dictionary.

    Repeated attributes such as 'tag' are stored as a semicolon-separated string.
    """

    attributes = {}

    for item in attribute_text.strip().split(";"):
        item = item.strip()

        if not item:
            continue

        parts = item.split(" ", 1)

        if len(parts) != 2:
            continue

        key, value = parts
        value = value.strip('"')

        if key in attributes:
            attributes[key] = f"{attributes[key]};{value}"
        else:
            attributes[key] = value

    return attributes


def load_gtf_exons(annotation_gtf, gene_name, transcript_id=None):
    """
    Load exon coordinates for one gene or transcript from an Ensembl/GENCODE GTF.

    Coordinates returned are 1-based, matching the GTF convention and the
    plotting functions used in this script.
    """

    annotation_gtf = Path(annotation_gtf)

    if not annotation_gtf.exists():
        raise FileNotFoundError(
            f"Annotation GTF file not found: {annotation_gtf}"
        )

    rows = []

    with open(annotation_gtf, "r", encoding="utf-8") as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attributes_text = fields

            if feature != "exon":
                continue

            attributes = parse_gtf_attributes(attributes_text)

            row_gene_name = attributes.get("gene_name")
            row_transcript_id = attributes.get("transcript_id")

            if gene_name is not None and row_gene_name != gene_name:
                continue

            if transcript_id is not None and row_transcript_id != transcript_id:
                continue

            rows.append({
                "chromosome": chrom,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "gene_name": row_gene_name,
                "gene_id": attributes.get("gene_id"),
                "transcript_id": row_transcript_id,
                "transcript_name": attributes.get("transcript_name"),
                "transcript_biotype": attributes.get("transcript_biotype"),
                "exon_number": attributes.get("exon_number"),
                "tag": attributes.get("tag"),
                "ccds_id": attributes.get("ccds_id"),
                "type": "exon",
            })

    exon_df = pd.DataFrame(rows)

    if exon_df.empty:
        raise ValueError(
            f"No exons found in {annotation_gtf} for gene_name={gene_name} "
            f"and transcript_id={transcript_id}."
        )

    return exon_df

def select_transcript_exons(exon_df, transcript_id=None):
    """
    Select one transcript from exon annotations.

    Priority:
    1. User-provided transcript_id.
    2. MANE_Select tag if present.
    3. Ensembl_canonical tag if present.
    4. CCDS-tagged transcript if present.
    5. Transcript with the highest number of exons.
    """

    selection_method = None

    if transcript_id is not None:
        selected_df = exon_df[
            exon_df["transcript_id"] == transcript_id
        ].copy()

        if selected_df.empty:
            raise ValueError(
                f"Transcript ID '{transcript_id}' was not found in the loaded exon annotations."
            )

        selection_method = "user_provided_transcript_id"

    else:
        candidate_df = exon_df.copy()

        mane_df = candidate_df[
            candidate_df["tag"].fillna("").str.contains("MANE_Select", na=False)
        ]

        if not mane_df.empty:
            candidate_df = mane_df
            selection_method = "MANE_Select"

        else:
            canonical_df = candidate_df[
                candidate_df["tag"].fillna("").str.contains("Ensembl_canonical", na=False)
            ]

            if not canonical_df.empty:
                candidate_df = canonical_df
                selection_method = "Ensembl_canonical"

            else:
                ccds_df = candidate_df[
                    candidate_df["tag"].fillna("").str.contains("CCDS", na=False)
                ]

                if not ccds_df.empty:
                    candidate_df = ccds_df
                    selection_method = "CCDS"

                else:
                    selection_method = "highest_exon_count"

        transcript_counts = (
            candidate_df
            .groupby("transcript_id")
            .size()
            .sort_values(ascending=False)
        )

        selected_transcript_id = transcript_counts.index[0]

        selected_df = exon_df[
            exon_df["transcript_id"] == selected_transcript_id
        ].copy()

    selected_df = selected_df.sort_values(
        by=["start", "end"]
    ).reset_index(drop=True)

    fallback_exon_numbers = pd.Series(
        selected_df.index + 1,
        index=selected_df.index
    ).astype(str)

    exon_numbers = selected_df["exon_number"].fillna(
        fallback_exon_numbers
    ).astype(str)

    selected_df["name"] = "exon " + exon_numbers

    first_row = selected_df.iloc[0]

    transcript_info = {
        "gene_name": first_row.get("gene_name"),
        "gene_id": first_row.get("gene_id"),
        "transcript_id": first_row.get("transcript_id"),
        "transcript_name": first_row.get("transcript_name"),
        "transcript_biotype": first_row.get("transcript_biotype"),
        "strand": first_row.get("strand"),
        "selection_method": selection_method,
        "n_exons": len(selected_df),
    }

    return selected_df, transcript_info

def build_exon_intron_annotation(selected_exons):
    """
    Build exon/intron annotation from one selected transcript.

    Exons come from the GTF. Introns are inferred as the regions between
    consecutive exons.
    """

    selected_exons = selected_exons.sort_values(
        by=["start", "end"]
    ).reset_index(drop=True)

    rows = []

    for i, exon in selected_exons.iterrows():
        rows.append({
            "name": exon["name"],
            "start": int(exon["start"]),
            "end": int(exon["end"]),
            "type": "exon",
        })

        if i < len(selected_exons) - 1:
            next_exon = selected_exons.iloc[i + 1]

            intron_start = int(exon["end"]) + 1
            intron_end = int(next_exon["start"]) - 1

            if intron_start <= intron_end:
                rows.append({
                    "name": f"intron {i + 1}-{i + 2}",
                    "start": intron_start,
                    "end": intron_end,
                    "type": "intron",
                })

    return pd.DataFrame(rows)

def resolve_annotation(args):
    """
    Resolve exon/intron annotation for plotting.

    Priority:
    1. If annotation_gtf and gene_name are provided, load annotation from GTF.
    2. Otherwise, fall back to the current SEC23B hardcoded annotation.
    """

    if args.annotation_gtf is not None:
        if args.gene_name is None:
            raise ValueError(
                "--annotation-gtf was provided, but --gene-name is missing. "
                "Provide a gene name such as SEC23B or DLG1."
            )

        exon_df = load_gtf_exons(
            annotation_gtf=args.annotation_gtf,
            gene_name=args.gene_name,
            transcript_id=args.transcript_id,
        )

        selected_exons, transcript_info = select_transcript_exons(
            exon_df=exon_df,
            transcript_id=args.transcript_id,
        )

        annotation_df = build_exon_intron_annotation(selected_exons)
        print(
            "Gene annotation loaded from GTF: "
            f"{transcript_info['gene_name']} | "
            f"{transcript_info['transcript_name']} | "
            f"{transcript_info['transcript_id']} | "
            f"selection: {transcript_info['selection_method']}"
        )   

        return annotation_df, transcript_info
    
    transcript_info = {
    "gene_name": None,
    "gene_id": None,
    "transcript_id": None,
    "transcript_name": None,
    "transcript_biotype": None,
    "strand": None,
    "selection_method": "no_annotation_provided",
    "n_exons": None,
    }

    print(
        "No gene annotation was provided. "
        "Figures will be generated without exon/intron annotation."
    )

    return None, transcript_info




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

def junctiondata_to_delta_dataframe(ref_junctions, alt_junctions, strand, region_start, region_end):
    """
    Convert AlphaGenome splice_junctions REF/ALT into a local dataframe
    with delta values.

    One row = one splice junction.
    """

    rows = []

    for junction_index, (junction, ref_values, alt_values) in enumerate(zip(
        ref_junctions.junctions,
        ref_junctions.values,
        alt_junctions.values,
    )):
        if junction.strand != strand:
            continue

        junction_start_1based = junction.start + 1
        junction_end_1based = junction.end

        overlaps_region = (
            junction_start_1based <= region_end
            and junction_end_1based >= region_start
        )

        if not overlaps_region:
            continue

        ref_value = float(ref_values[0])
        alt_value = float(alt_values[0])
        delta_alt_ref = alt_value - ref_value

        rows.append({
            "junction_index": junction_index,
            "chromosome": junction.chromosome,
            "junction_start_1based": junction_start_1based,
            "junction_end_1based": junction_end_1based,
            "strand": junction.strand,
            "ref_value": ref_value,
            "alt_value": alt_value,
            "delta_alt_ref": delta_alt_ref,
            "abs_delta": abs(delta_alt_ref),
        })

    df = pd.DataFrame(rows)

    if df.empty:
        return df

    return df.sort_values("abs_delta", ascending=False).reset_index(drop=True)

def select_junctions_by_directional_delta(
    junction_df,
    min_positive_delta,
    max_negative_delta,
):
    """
    Select junctions by directional ALT-REF delta.

    Increased junctions:
        delta_alt_ref >= min_positive_delta

    Decreased junctions:
        delta_alt_ref <= max_negative_delta
    """

    if junction_df.empty:
        return junction_df

    selected_df = junction_df[
        (junction_df["delta_alt_ref"] >= min_positive_delta)
        |
        (junction_df["delta_alt_ref"] <= max_negative_delta)
    ].copy()

    if selected_df.empty:
        return selected_df

    selected_df["change_type"] = "unchanged"

    selected_df.loc[
        selected_df["delta_alt_ref"] >= min_positive_delta,
        "change_type"
    ] = "increased"

    selected_df.loc[
        selected_df["delta_alt_ref"] <= max_negative_delta,
        "change_type"
    ] = "decreased"

    selected_df = selected_df.sort_values(
        by="abs_delta",
        ascending=False
    ).reset_index(drop=True)

    return selected_df


def plot_junction_arc(
    ax,
    start,
    end,
    y_base,
    height,
    color,
    linewidth,
    label_text=None,
):
    """
    Draw a splice-junction arc between two genomic coordinates.
    """

    x_mid = (start + end) / 2
    width = end - start

    arc = Arc(
        (x_mid, y_base),
        width=width,
        height=height,
        angle=0,
        theta1=0,
        theta2=180,
        color=color,
        linewidth=linewidth,
    )

    ax.add_patch(arc)

    if label_text is not None:
        ax.text(
            x_mid,
            y_base + height / 2 + 0.12,
            label_text,
            ha="center",
            va="bottom",
            fontsize=8,
            color=color,
            fontweight="bold",
        )


def plot_selected_splice_junctions(
    selected_junctions,
    region_start,
    region_end,
    variant_pos,
    output_dir,
    title,
    strand,
    annotation_df=None,
    transcript_info=None
):
    """
    Plot selected splice junctions using directional ALT-REF filtering.

    Decreased junctions are shown as REF-associated arcs.
    Increased junctions are shown as ALT-associated arcs.
    Arc thickness is proportional to abs(delta_alt_ref).
    Arc height is proportional to junction length.
    """

    if selected_junctions.empty:
        print("No splice junctions passed the selected thresholds.")
        return None

    fig, (ax_annotation, ax_junctions) = plt.subplots(
        nrows=2,
        figsize=(12, 4.8),
        sharex=True,
        gridspec_kw={"height_ratios": [1, 4]}
    )

    ax_junctions.set_xlim(region_start, region_end)
    ax_junctions.set_ylim(-0.05, 1.35)

    ax_junctions.set_title(f"{title} — selected {SPLICE_JUNCTIONS_LABEL}\n{strand} strand")

    ax_junctions.set_xticks([])
    ax_junctions.set_yticks([])
    ax_junctions.set_xlabel("")
    ax_junctions.set_ylabel("")

    for spine in ax_junctions.spines.values():
        spine.set_visible(False)

    max_abs_delta = selected_junctions["abs_delta"].max()
    if max_abs_delta == 0:
        max_abs_delta = 1

    max_length = (
        selected_junctions["junction_end_1based"]
        - selected_junctions["junction_start_1based"]
    ).max()
    if max_length == 0:
        max_length = 1


    for _, row in selected_junctions.iterrows():
        start = row["junction_start_1based"]
        end = row["junction_end_1based"]
        junction_length = end - start

        arc_height = 0.35 + 0.65 * (junction_length / max_length)
        linewidth = 1.5 + 4.5 * (row["abs_delta"] / max_abs_delta)
        delta_label = f"Δ {row['delta_alt_ref']:+.2f}"
     
       

        if row["change_type"] == "decreased":
            plot_junction_arc(
                ax=ax_junctions,
                start=start,
                end=end,
                y_base=0,
                height=arc_height,
                color=REF_COLOR,
                linewidth=linewidth,
                label_text=delta_label,
            )

        elif row["change_type"] == "increased":
            plot_junction_arc(
                ax=ax_junctions,
                start=start,
                end=end,
                y_base=0,
                height=arc_height,
                color=ALT_COLOR,
                linewidth=linewidth,
                label_text=delta_label,
            )

    if variant_pos is not None:
        ax_junctions.axvline(
            variant_pos,
            color=VARIANT_COLOR,
            linestyle="-.",
            linewidth=1.2,
        )

        ax_junctions.text(
            variant_pos - 150,
            1.20,
            f"KI\n({variant_pos:,})",
            color=VARIANT_COLOR,
            ha="right",
            va="top",
            fontsize=9,
            fontweight="bold",
        )

    legend_elements = [
        Line2D([0], [0], color=REF_COLOR, linewidth=2.5, label="Decreased in ALT"),
        Line2D([0], [0], color=ALT_COLOR, linewidth=2.5, label="Increased in ALT"),
    ]

    ax_junctions.legend(
        handles=legend_elements,
        loc="upper right",
        frameon=True,
    )

    if annotation_df is not None:
        plot_mane_annotation_track(
            ax=ax_annotation,
            annotation_df=annotation_df,
            region_start=region_start,
            region_end=region_end,
            variant_pos=variant_pos,
            transcript_info=transcript_info,
        )
    else:
        ax_annotation.axis("off")

    fig.tight_layout()

    strand_suffix = "plus" if strand == "+" else "minus"
    output_path = output_dir / f"04_splice_junctions_{strand_suffix}.png"

    fig.savefig(output_path, dpi=300)
    plt.close(fig)

    print(f"{SPLICE_JUNCTIONS_LABEL} figure saved to: {output_path}")

    return output_path

def write_visualization_runlog(
    args,
    output_dir,
    prediction_pkl,
    input_runlog_path,
    region_start,
    region_end,
    region_source,
    variant_pos,
    transcript_info,
    selected_junctions,
    generated_figure_paths,
):
    """
    Write runlog for visualize-prediction.py execution.
    """

    runlog_path = output_dir / "visualization_runlog.txt"

    with open(runlog_path, "w", encoding="utf-8") as runlog:
        runlog.write("Experiment exp03 — AlphaGenome Prediction Visualization\n")
        runlog.write("=" * 60 + "\n\n")

        runlog.write(f"Execution date: {datetime.now()}\n\n")

        runlog.write("Input files\n")
        runlog.write("-" * 20 + "\n")
        runlog.write(f"prediction_pkl: {prediction_pkl}\n")
        runlog.write(f"input_runlog: {input_runlog_path}\n\n")
        runlog.write("Gene annotation\n")
        runlog.write("-" * 20 + "\n")

        for key, value in transcript_info.items():
            runlog.write(f"{key}: {value}\n")

        runlog.write(f"annotation_gtf: {args.annotation_gtf}\n")
        runlog.write(f"gene_name_argument: {args.gene_name}\n")
        runlog.write(f"transcript_id_argument: {args.transcript_id}\n")
        runlog.write("\n")
        
        runlog.write("Visualization parameters\n")
        runlog.write("-" * 20 + "\n")
        runlog.write(f"output_dir: {output_dir}\n")
        runlog.write(f"title: {args.title}\n")
        runlog.write(f"strand: {args.strand}\n")
        runlog.write(f"region_start_1based: {region_start}\n")
        runlog.write(f"region_end_1based: {region_end}\n")
        runlog.write(f"region_source: {region_source}\n")
        runlog.write(f"variant_pos_1based: {variant_pos}\n")
        runlog.write(f"view_flank: {args.view_flank}\n\n")

        runlog.write("Splice junction selection parameters\n")
        runlog.write("-" * 20 + "\n")
        runlog.write(f"junction_min_positive_delta: {args.junction_min_positive_delta}\n")
        runlog.write(f"junction_max_negative_delta: {args.junction_max_negative_delta}\n")
        runlog.write(f"selected_junction_count: {len(selected_junctions)}\n\n")

        runlog.write("Selected splice junctions\n")
        runlog.write("-" * 20 + "\n")
        if selected_junctions.empty:
            runlog.write("No splice junctions selected.\n")
        else:
            runlog.write(selected_junctions.to_string(index=False))
            runlog.write("\n")
        runlog.write("\n")

        runlog.write("Generated figures\n")
        runlog.write("-" * 20 + "\n")
        for path in generated_figure_paths:
            if path is not None:
                runlog.write(f"{path.name}\n")

    print(f"Visualization runlog saved to: {runlog_path}")

    return runlog_path

def main():
    args = parse_arguments()

    prediction_pkl = validate_existing_file(
    args.prediction_pkl,
    "--prediction-pkl"
    )

    prediction_dir = prediction_pkl.parent
    
    runlog_path = resolve_runlog_path(
    args,
    prediction_dir
    )
    validate_annotation_arguments(args)
    
    output_dir = Path(args.output_dir) if args.output_dir else prediction_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)

    generated_figure_paths = []

    prediction = load_prediction(prediction_pkl)
    available_outputs = get_available_prediction_outputs(prediction)
    
    region_start, region_end, region_source = resolve_region(
    args,
    runlog_path
    )

    variant_pos = resolve_variant_position(args, runlog_path)
    
    validate_variant_inside_region(
    variant_pos,
    region_start,
    region_end
    )

    print("Prediction loaded successfully.")
    print(f"Prediction file: {prediction_pkl}")
    print(f"Output directory: {output_dir}")
    print(f"Runlog: {runlog_path}")
    print(f"Region: {region_start}-{region_end}")
    print(f"Region source: {region_source}")
    print(f"Variant position: {variant_pos}")
    print(f"Strand: {args.strand}")

    print("\nAvailable outputs required by the visualization script:")

    for output_name, is_available in available_outputs.items():
        status = "available" if is_available else "not available"
        print(f"- {output_name}: {status}")
    annotation_df, transcript_info = resolve_annotation(args)

    if args.strand in ["+", "-"]:
        splice_site_indices = get_splice_site_track_indices(args.strand)
        print("\nSplice site track indices:")
        print(splice_site_indices)

    if available_outputs["rna_seq"]:    
        rna_seq_path = plot_rna_seq_ref_alt_and_delta(
            prediction=prediction,
            region_start=region_start,
            region_end=region_end,
            variant_pos=variant_pos,
            output_dir=output_dir,
            title=args.title,
            annotation_df=annotation_df,
            transcript_info=transcript_info
        )
        generated_figure_paths.append(rna_seq_path)
    else:
        print("Skipping RNA-seq figure: rna_seq output is not available in prediction.pkl.")

    if available_outputs["splice_sites"]:
        if args.strand in ["+", "-"] and available_outputs["splice_junctions"]:
            splice_sites_path = plot_splice_sites(
                prediction=prediction,
                region_start=region_start,
                region_end=region_end,
                variant_pos=variant_pos,
                output_dir=output_dir,
                title=args.title,
                strand=args.strand,
                annotation_df=annotation_df,
                transcript_info=transcript_info
            )
            generated_figure_paths.append(splice_sites_path)

        else:
            for strand in ["+", "-"]:
                splice_sites_path = plot_splice_sites(
                    prediction=prediction,
                    region_start=region_start,
                    region_end=region_end,
                    variant_pos=variant_pos,
                    output_dir=output_dir,
                    title=args.title,
                    strand=strand,
                    annotation_df=annotation_df,
                    transcript_info=transcript_info
                )
                generated_figure_paths.append(splice_sites_path)
    else:
        print("Skipping splice_sites figure: splice_sites output is not available in prediction.pkl.")
    
    if available_outputs["splice_site_usage"]:
        if args.strand in ["+", "-"]:
            splice_site_usage_path = plot_splice_site_usage(
                prediction=prediction,
                region_start=region_start,
                region_end=region_end,
                variant_pos=variant_pos,
                output_dir=output_dir,
                title=args.title,
                strand=args.strand,
                annotation_df=annotation_df,
                transcript_info=transcript_info
            )
            generated_figure_paths.append(splice_site_usage_path)

        else:
            for strand in ["+", "-"]:
                splice_site_usage_path = plot_splice_site_usage(
                    prediction=prediction,
                    region_start=region_start,
                    region_end=region_end,
                    variant_pos=variant_pos,
                    output_dir=output_dir,
                    title=args.title,
                    strand=strand,
                    annotation_df=annotation_df,
                    transcript_info=transcript_info
                )
                generated_figure_paths.append(splice_site_usage_path)
    else:
        print("Skipping splice_site_usage figure: splice_site_usage output is not available in prediction.pkl.")
    
    selected_junctions = pd.DataFrame()

    if args.strand in ["+", "-"]:
        junction_df = junctiondata_to_delta_dataframe(
            ref_junctions=prediction.reference.splice_junctions,
            alt_junctions=prediction.alternate.splice_junctions,
            strand=args.strand,
            region_start=region_start,
            region_end=region_end,
        )

        selected_junctions = select_junctions_by_directional_delta(
            junction_df=junction_df,
            min_positive_delta=args.junction_min_positive_delta,
            max_negative_delta=args.junction_max_negative_delta,
        )

        print("\nSelected splice junctions by directional delta:")

        if selected_junctions.empty:
            print(
                f"No splice junctions selected with "
                f"delta >= {args.junction_min_positive_delta} "
                f"or delta <= {args.junction_max_negative_delta}"
            )
        else:
            print(selected_junctions.to_string(index=False))

        splice_junctions_path = plot_selected_splice_junctions(
            selected_junctions=selected_junctions,
            region_start=region_start,
            region_end=region_end,
            variant_pos=variant_pos,
            output_dir=output_dir,
            title=args.title,
            strand=args.strand,
            annotation_df=annotation_df,
            transcript_info=transcript_info
        )
        generated_figure_paths.append(splice_junctions_path)

    else:
        if args.strand not in ["+", "-"]:
            print("\nSkipping splice_junctions plot: use --strand + or --strand -.")
        elif not available_outputs["splice_junctions"]:
            print("\nSkipping splice_junctions plot: splice_junctions output is not available in prediction.pkl.")
        

    write_visualization_runlog(
        args=args,
        output_dir=output_dir,
        prediction_pkl=prediction_pkl,
        input_runlog_path=runlog_path,
        region_start=region_start,
        region_end=region_end,
        region_source=region_source,
        variant_pos=variant_pos,
        transcript_info=transcript_info,
        selected_junctions=selected_junctions,
        generated_figure_paths=generated_figure_paths,
    )   
            
if __name__ == "__main__":
    main()