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

fig, ax = plt.subplots(figsize=(10, 3.8))


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
    VARIANT_POS + 140,
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

ax.set_xlim(197076044, 197086544)
ax.ticklabel_format(style="plain", axis="x")
ax.get_xaxis().get_major_formatter().set_useOffset(False)
ax.set_xlabel("Genomic position (GRCh38)")
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

plt.savefig(
    out_png,
    dpi=300,
    bbox_inches="tight",
)
plt.show()

