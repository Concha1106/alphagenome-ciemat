#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 19:24:41 2026
Description: Script to graph the rna_seq parameter of figure 3b

@author: concha
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

REF_COLOR = "#1f77b4"
ALT_COLOR = "#7a0019"
VARIANT_COLOR = "#D4A017"

TITLE = "DLG1 — Artery tibial:\npredicted exon skipping"

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

        
        
# 1) Paths

project_root = Path(__file__).resolve().parents[2]
results_dir = project_root / "results/2026-04-10-exp02-exon-skipping"

rna_file = results_dir / "dlg1_rna_seq_ref_vs_alt.tsv"
out_png = results_dir /"dlg1_fig3b_rnaseq_improved.png"

# 2) Load data

df = pd.read_csv(rna_file, sep="\t")

print(df.columns)
print(df.head())

# 3) Smooth signal

window = 150  # probar con diferntes valores, ver con Andrés (rolling calcula la media de los puntos vecinos)

df["ref_smooth"] = df["ref_value"].rolling(window=window, center=True).mean()
df["alt_smooth"] = df["alt_value"].rolling(window=window, center=True).mean()


# 3) Plot

fig, (ax, ax_annot) = plt.subplots(
    nrows=2,
    ncols=1,
    figsize=(10, 4.8),
    sharex=True,
    gridspec_kw={"height_ratios": [3.8, 0.7]},
)

ax.plot(
    df["position"],
    df["ref_smooth"],
    label="REF",
    color=REF_COLOR,
    linewidth=2,
)

ax.plot(
    df["position"],
    df["alt_smooth"],
    label="ALT",
    color=ALT_COLOR,
    linewidth=2,
)

ax.axvline(
    VARIANT_POS,
    color=VARIANT_COLOR,
    linewidth=1.2,
    linestyle="--",
)

ax.text(
    VARIANT_POS - 350,
   0.96,
    "4 bp deletion\n197081044",
    fontsize=8,
    color=VARIANT_COLOR,
    ha="left",
    va="top",
    fontweight="bold",
)

# 4) Formatting

ax.set_xlim(REGION_END, REGION_START)
ax.ticklabel_format(style="plain", axis="x")
ax.get_xaxis().get_major_formatter().set_useOffset(False)

ax.set_ylabel("Predicted\nRNA-seq signal")
ax.set_xlabel("")

ax.set_title(
    TITLE,
    fontsize=11,
    fontweight="bold",
)

ax.legend(
    frameon=False,
    loc="center left",
    bbox_to_anchor=(1.01, 0.92),
)


ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

def plot_dlg1_annotation(ax):
    ax.set_xlim(REGION_END, REGION_START)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylabel("DLG1", rotation=0, ha="right", va="center")

    # intron line
    ax.hlines(
        y=0.5,
        xmin=min(exon["start"] for exon in DLG1_EXONS),
        xmax=max(exon["end"] for exon in DLG1_EXONS),
        color=INTRON_COLOR,
        linewidth=1.2,
    )

    # exons
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
        
plot_dlg1_annotation(ax_annot)

plt.tight_layout()
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.show()
