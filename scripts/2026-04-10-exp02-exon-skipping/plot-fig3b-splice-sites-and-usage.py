#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 01 2026
Description: Script to graph the splice site and splice sites usage parameters of figure 3b

@author: concha
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

REF_COLOR = "#1f77b4"
ALT_COLOR = "#7a0019"
VARIANT_COLOR = "#D4A017"

TITLE = "DLG1 | Artery tibial\nPredicted splice sites and usage"

VARIANT_POS = 197081044

REGION_START = 197076044
REGION_END = 197086544

EXON_COLOR = "#1f77b4"
INTRON_COLOR = "#555555"

DLG1_EXONS = [
    {"name": "E16", "start": 197085580, "end": 197085756},
    {"name": "E17", "start": 197081051, "end": 197081117},
    {"name": "E18", "start": 197076586, "end": 197076685},
]

def plot_dlg1_annotation(ax):
    ax.set_xlim(REGION_END, REGION_START)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylabel("DLG1", rotation=0, ha="right", va="center")

    ax.hlines(
        y=0.5,
        xmin=min(exon["start"] for exon in DLG1_EXONS),
        xmax=max(exon["end"] for exon in DLG1_EXONS),
        color=INTRON_COLOR,
        linewidth=1.2,
    )

    for exon in DLG1_EXONS:
        width = exon["end"] - exon["start"] + 1

        ax.broken_barh(
            [(exon["start"], width)],
            (0.42, 0.22),
            facecolors=EXON_COLOR,
            edgecolors=EXON_COLOR,
        )

        ax.text(
            exon["start"] + width / 2,
            0.72,
            exon["name"],
            ha="center",
            va="bottom",
            fontsize=8,
            fontweight="bold",
        )

    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.set_xlabel("Genomic position, chr3:197086544-197076044 (GRCh38)")

# 1) Paths
project_root = Path("~/proyectos_UAX/alphagenome-ciemat").expanduser()
results_dir  = project_root / "results/2026-04-10-exp02-exon-skipping"

sites_file = results_dir / "dlg1_splice_sites_ref_vs_alt.tsv"
usage_file = results_dir / "dlg1_splice_site_usage_ref_vs_alt.tsv"
out_png    = results_dir / "dlg1_fig3b_sites_usage_improved.png"

# 2) Load data
df_sites = pd.read_csv(sites_file, sep="\t")
df_usage = pd.read_csv(usage_file, sep="\t")

# 3) Smooth signal (mismo criterio que RNA-seq)
window = 150
for df in [df_sites, df_usage]:
    df["ref_smooth"] = df["ref_value"].rolling(window=window, center=True).mean()
    df["alt_smooth"] = df["alt_value"].rolling(window=window, center=True).mean()

# Splice sites tiene donor y acceptor — los separamos
donor    = df_sites[df_sites["track"] == "donor"]
acceptor = df_sites[df_sites["track"] == "acceptor"]

fig, (ax_sites, ax_usage, ax_annot) = plt.subplots(
    nrows=3,
    ncols=1,
    figsize=(10, 6.0),
    sharex=True,
    gridspec_kw={"height_ratios": [2.4, 2.4, 0.7]},
)

# --- Splice sites ---
ax_sites.plot(donor["position"], donor["ref_smooth"], color=REF_COLOR, linewidth=1.5)
ax_sites.plot(donor["position"], donor["alt_smooth"], color=ALT_COLOR, linewidth=1.5)
ax_sites.plot(acceptor["position"], acceptor["ref_smooth"], color=REF_COLOR, linewidth=1.5)
ax_sites.plot(acceptor["position"], acceptor["alt_smooth"], color=ALT_COLOR, linewidth=1.5)

ax_sites.axvline(VARIANT_POS, color=VARIANT_COLOR, linewidth=1.2, linestyle="--")

ax_sites.set_title(
    TITLE,
    fontsize=11,
    fontweight="bold",
)
ax_sites.set_ylabel("Splice sites")
ax_sites.set_yticks([])

# --- Splice site usage ---
ax_usage.plot(df_usage["position"], df_usage["ref_smooth"], color=REF_COLOR, linewidth=1.5, label="REF")
ax_usage.plot(df_usage["position"], df_usage["alt_smooth"], color=ALT_COLOR, linewidth=1.5, label="ALT")

ax_usage.axvline(VARIANT_POS, color=VARIANT_COLOR, linewidth=1.2, linestyle="--")
ax_usage.set_xlabel("")
ax_usage.set_yticks([])

ax_sites.text(
    VARIANT_POS - 350,
    ax_sites.get_ylim()[1]*0.92,
    "4 bp deletion\n197081044",
    fontsize=8,
    color=VARIANT_COLOR,
    ha="left",
    va="top",
    fontweight="bold",
)

# 5) Formatting
ax_usage.set_xlim(REGION_END, REGION_START)
ax_usage.ticklabel_format(style="plain", axis="x")
ax_usage.get_xaxis().get_major_formatter().set_useOffset(False)
ax_usage.set_xlabel("Genomic position (GRCh38)")

leyenda = [
    Line2D([0], [0], color=REF_COLOR, linewidth=2.5, label="REF"),
    Line2D([0], [0], color=ALT_COLOR, linewidth=2.5, label="ALT"),
]

ax_usage.legend(
    handles=leyenda,
    frameon=False,
    loc="center left",
    bbox_to_anchor=(1.01, 1.00),
)

for ax in [ax_sites, ax_usage]:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

plot_dlg1_annotation(ax_annot)
plt.tight_layout()
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.show()