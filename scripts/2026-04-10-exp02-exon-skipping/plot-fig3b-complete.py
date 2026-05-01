#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 12:42:56 2026
Script to replicate the complete graph of figure 3b
@author: concha
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from matplotlib.lines import Line2D

# 1) Paths
project_root   = Path("~/Desktop/test-zip-recovery/alphagenome-ciemat").expanduser()
results_dir    = project_root / "results/2026-04-10-exp02-exon-skipping"

junctions_file = results_dir / "dlg1_splice_junctions_ref_vs_alt.tsv"
rna_file       = results_dir / "dlg1_rna_seq_ref_vs_alt.tsv"
sites_file     = results_dir / "dlg1_splice_sites_ref_vs_alt.tsv"
usage_file     = results_dir / "dlg1_splice_site_usage_ref_vs_alt.tsv"
out_png        = results_dir / "dlg1_fig3b_complete.png"

# 2) Load data
df_junc  = pd.read_csv(junctions_file, sep="\t")
df_rna   = pd.read_csv(rna_file,       sep="\t")
df_sites = pd.read_csv(sites_file,     sep="\t")
df_usage = pd.read_csv(usage_file,     sep="\t")

# 3) Smooth RNA-seq and sites/usage signals
window = 150
for df in [df_rna, df_sites, df_usage]:
    df["ref_smooth"] = df["ref_value"].rolling(window=window, center=True).mean()
    df["alt_smooth"] = df["alt_value"].rolling(window=window, center=True).mean()

donor    = df_sites[df_sites["track"] == "donor"]
acceptor = df_sites[df_sites["track"] == "acceptor"]

# 4) Select junctions of interest
junctions_of_interest = pd.DataFrame({
    "start": [197076685, 197081117, 197076685],
    "end":   [197081050, 197085579, 197085579],
    "panel": ["REF",     "REF",     "ALT"],
})
plot_junc = junctions_of_interest.merge(df_junc, on=["start", "end"], how="left")

# 5) Figure with 4 panels
# height_ratios: junctions panel taller, the rest equal
fig, (ax_junc, ax_sites, ax_usage, ax_rna) = plt.subplots(
    nrows=4, ncols=1,
    figsize=(10, 9),
    sharex=True,
    gridspec_kw={"height_ratios": [3, 1, 1, 1.5], "hspace": 0.08},
)

fig.suptitle(
    "DLG1 — chr3:197081044 TACTC>T  |  Tissue: artery (tibial)",
    fontsize=11, fontweight="bold"
)

VARIANT_POS = 197081044

# 6) Panel 1: Splice junctions (arcs)
panel_y = {"REF": 1.00, "ALT": 0.35}

for _, row in plot_junc.iterrows():
    x_mid  = (row["start"] + row["end"]) / 2
    width  = row["end"] - row["start"]
    y      = panel_y[row["panel"]]
    height = 0.4 if row["panel"] == "REF" else 0.55
    color  = "cornflowerblue" if row["panel"] == "REF" else "firebrick"
    lw     = 1 + abs(row["delta_alt_ref"]) * 4

    arc = Arc(
        (x_mid, y),
        width=width, height=height,
        theta1=0, theta2=180,
        color=color, linewidth=lw,
    )
    ax_junc.add_patch(arc)

    ax_junc.text(
        x_mid, y + height / 2 + 0.08,
        f"{abs(row['delta_alt_ref']):.2f}",
        ha="center", va="bottom", fontsize=10,
    )

ax_junc.axvline(VARIANT_POS, color="gold", linewidth=1.5, linestyle="--")
ax_junc.set_ylim(0.2, 2.0)
ax_junc.set_yticks([])
ax_junc.set_ylabel("Splice\njunctions")

leyenda = [
    Line2D([0], [0], color="cornflowerblue", linewidth=2.5, label="REF"),
    Line2D([0], [0], color="firebrick",      linewidth=2.5, label="ALT"),
]
ax_junc.legend(handles=leyenda, loc="upper right", frameon=False)

# 7) Panel 2: Splice sites
ax_sites.plot(donor["position"],    donor["ref_smooth"],    color="cornflowerblue", linewidth=2)
ax_sites.plot(donor["position"],    donor["alt_smooth"],    color="firebrick",      linewidth=2)
ax_sites.plot(acceptor["position"], acceptor["ref_smooth"], color="cornflowerblue", linewidth=2)
ax_sites.plot(acceptor["position"], acceptor["alt_smooth"], color="firebrick",      linewidth=2)
ax_sites.axvline(VARIANT_POS, color="gold", linewidth=1.5, linestyle="--")
ax_sites.set_ylabel("Splice\nsites")
ax_sites.set_yticks([])

# 8) Panel 3: Splice site usage
ax_usage.plot(df_usage["position"], df_usage["ref_smooth"], color="cornflowerblue", linewidth=2)
ax_usage.plot(df_usage["position"], df_usage["alt_smooth"], color="firebrick",      linewidth=2)
ax_usage.axvline(VARIANT_POS, color="gold", linewidth=1.5, linestyle="--")
ax_usage.set_ylabel("Splice site\nusage")
ax_usage.set_yticks([])

# 9) Panel 4: RNA-seq
ax_rna.plot(df_rna["position"], df_rna["ref_smooth"], color="cornflowerblue", linewidth=2)
ax_rna.plot(df_rna["position"], df_rna["alt_smooth"], color="firebrick",      linewidth=2)
ax_rna.axvline(VARIANT_POS, color="gold", linewidth=1.5, linestyle="--")
ax_rna.set_ylabel("RNA-seq\nsignal")
ax_rna.set_yticks([])

# 10) Formatting (only on bottom panel since sharex=True)
ax_rna.set_xlim(197076044, 197086544)
ax_rna.ticklabel_format(style="plain", axis="x")
ax_rna.get_xaxis().get_major_formatter().set_useOffset(False)
ax_rna.set_xlabel("Genomic position (chr3, hg38)")

plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.show()
