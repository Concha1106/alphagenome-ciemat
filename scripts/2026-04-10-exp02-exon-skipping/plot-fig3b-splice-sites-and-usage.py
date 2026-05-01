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

# 1) Paths
project_root = Path("~/Desktop/test-zip-recovery/alphagenome-ciemat").expanduser()
results_dir  = project_root / "results/2026-04-10-exp02-exon-skipping"

sites_file = results_dir / "dlg1_splice_sites_ref_vs_alt.tsv"
usage_file = results_dir / "dlg1_splice_site_usage_ref_vs_alt.tsv"
out_png    = results_dir / "dlg1_fig3b_sites_usage.png"

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

# 4) Plot
fig, (ax_sites, ax_usage) = plt.subplots(
    nrows=2, ncols=1,
    figsize=(10, 5),
    sharex=True,
)

# --- Splice sites ---## hacer con la media de donor y acceptor?
ax_sites.plot(donor["position"],    donor["ref_smooth"],    color="cornflowerblue", linewidth=2)
ax_sites.plot(donor["position"],    donor["alt_smooth"],    color="firebrick",      linewidth=2)
ax_sites.plot(acceptor["position"], acceptor["ref_smooth"], color="cornflowerblue", linewidth=2)
ax_sites.plot(acceptor["position"], acceptor["alt_smooth"], color="firebrick",      linewidth=2)

ax_sites.axvline(197081044, color="gold", linewidth=1.5, linestyle="--")
ax_sites.set_ylabel("Splice sites")
ax_sites.set_yticks([])

# --- Splice site usage ---
ax_usage.plot(df_usage["position"], df_usage["ref_smooth"], color="cornflowerblue", linewidth=2, label="REF")
ax_usage.plot(df_usage["position"], df_usage["alt_smooth"], color="firebrick",      linewidth=2, label="ALT")

ax_usage.axvline(197081044, color="gold", linewidth=1.5, linestyle="--")
ax_usage.set_ylabel("Splice site usage")
ax_usage.set_yticks([])

# 5) Formatting
ax_usage.set_xlim(197076044, 197086544)
ax_usage.ticklabel_format(style="plain", axis="x")
ax_usage.get_xaxis().get_major_formatter().set_useOffset(False)
ax_usage.set_xlabel("Genomic position")

leyenda = [
    Line2D([0], [0], color="cornflowerblue", linewidth=2.5, label="REF"),
    Line2D([0], [0], color="firebrick",      linewidth=2.5, label="ALT"),
]
ax_usage.legend(handles=leyenda, loc="upper right", frameon=False)

plt.tight_layout()
plt.savefig(out_png, dpi=300)
plt.show()