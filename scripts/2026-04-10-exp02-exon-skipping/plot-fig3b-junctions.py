#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 14:29:38 2026

@author: e6260
"""

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Arc


# 1) Paths

project_root = Path("~/Desktop/test-zip-recovery/alphagenome-ciemat").expanduser()
results_dir = project_root / "results/2026-04-10-exp02-exon-skipping"

junctions_file = results_dir / "dlg1_splice_junctions_ref_vs_alt.tsv"

out_png = results_dir / "dlg1_fig3b_arcs_only.png"


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

fig, ax = plt.subplots(figsize=(10, 3))

panel_y = {
    "REF": 0.8,
    "ALT": 1.8,
}

for _, row in plot_junc.iterrows():
    x_mid = (row["start"] + row["end"]) / 2
    width = row["end"] - row["start"]
    y = panel_y[row["panel"]]

    arc = Arc(
        (x_mid, y),
        width=width,
        height=0.8,
        theta1=0,
        theta2=180,
        linewidth=2,
    )

    ax.add_patch(arc)


# 6) Basic formatting

ax.set_xlim(197075000, 197087000)
ax.set_ylim(0.2, 2.6)
ax.set_yticks([])

plt.tight_layout()
plt.savefig(out_png, dpi=300)
plt.show()

