#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 14:29:38 2026
Description: Script to graph the splices junction parameter of figure 3b

@author: e6260
"""

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from matplotlib.lines import Line2D

REF_COLOR = "#1f77b4"
ALT_COLOR = "#7a0019"
VARIANT_COLOR = "#D4A017"

TITLE = "DLG1 | Artery tibial\nPredicted splice junctions"

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
results_dir = project_root / "results/2026-04-10-exp02-exon-skipping"

junctions_file = results_dir / "dlg1_splice_junctions_ref_vs_alt.tsv"

out_png = results_dir / "dlg1_fig3b_arcs_improved.png"


# 2) Load junction table

df_junc = pd.read_csv(junctions_file, sep="\t")


# 3) Select junctions of interest

junctions_of_interest = pd.DataFrame({
    "start": [197076685, 197081117, 197076685],
    "end":   [197081050, 197085579, 197085579],
    "panel": ["REF", "REF", "ALT"],
})


# 4) Add predicted values from AlphaGenome output table

plot_junc = junctions_of_interest.merge(
    df_junc,
    on=["start", "end"],
    how="left",
)


# 5) Draw only arcs

fig, (ax, ax_annot) = plt.subplots(
    nrows=2,
    ncols=1,
    figsize=(10, 4.8),
    sharex=True,
    gridspec_kw={"height_ratios": [3.8, 0.7]},
)


panel_y = {
    "REF": 1.00,
    "ALT": 0.35,
}

for _, row in plot_junc.iterrows():
    x_mid = (row["start"] + row["end"]) / 2
    width = row["end"] - row["start"]
    y = panel_y[row["panel"]]

    height = 0.4 if row["panel"] == "REF" else 0.55
    color = REF_COLOR if row["panel"] == "REF" else ALT_COLOR
    lw = 1.2 + abs(row["delta_alt_ref"])*3

    
    arc = Arc(
        (x_mid, y),
        width=width,
        height=height,
        theta1=0,
        theta2=180,
        color=color,
        linewidth=lw,
    )

    ax.add_patch(arc)
   
    ax.text(
    x_mid,
    y + height / 2 + 0.08,
    f"{abs(row['delta_alt_ref']):.2f}",
    ha="center",
    va="bottom",
    fontsize=9,
    fontweight="bold",
    )

ax.axvline(
    VARIANT_POS,
    color=VARIANT_COLOR,
    linewidth=1.2,
    linestyle="--",
)

ax.text(
    VARIANT_POS - 350,
    1.92,
    "4 bp deletion\n197081044",
    fontsize=8,
    color=VARIANT_COLOR,
    ha="left",
    va="top",
    fontweight="bold",
)


# 6) Basic formatting

ax.set_title(
    TITLE,
    fontsize=11,
    fontweight="bold",
)

ax.set_xlim(REGION_END, REGION_START)
ax.ticklabel_format(style="plain", axis="x")
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.set_xlabel("")
ax.set_ylim(0.2, 2.25)
ax.set_yticks([])

leyenda = [
    Line2D([0], [0], color=REF_COLOR, linewidth=2.5, label="REF"),
    Line2D([0], [0], color=ALT_COLOR, linewidth=2.5, label="ALT"),
]

ax.legend(
    handles=leyenda,
    frameon=False,
    loc="center left",
    bbox_to_anchor=(1.01,0.90),
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plot_dlg1_annotation(ax_annot)
plt.tight_layout()


plt.savefig(
    out_png,
    dpi=300,
    bbox_inches="tight",
)
plt.show()

