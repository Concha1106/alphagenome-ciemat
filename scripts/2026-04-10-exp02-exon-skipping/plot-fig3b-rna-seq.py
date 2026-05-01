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

# 1) Paths

project_root = Path("~/Desktop/test-zip-recovery/alphagenome-ciemat").expanduser()
results_dir = project_root / "results/2026-04-10-exp02-exon-skipping"

rna_file = results_dir / "dlg1_rna_seq_ref_vs_alt.tsv"
out_png = results_dir / "dlg1_fig3b_rnaseq_only.png"

# 2) Load data

df = pd.read_csv(rna_file, sep="\t")

print(df.columns)
print(df.head())

# 3) Smooth signal

window = 150  # probar con diferntes valores, ver con Andrés (rolling calcula la media de los puntos vecinos)

df["ref_smooth"] = df["ref_value"].rolling(window=window, center=True).mean()
df["alt_smooth"] = df["alt_value"].rolling(window=window, center=True).mean()


# 3) Plot  ##pendiente ver que dejo

fig, ax = plt.subplots(figsize=(10, 3.5))  # por qué estos valores? cambiar aqui?

ax.plot(df["position"], df["ref_smooth"], label="REF", color="cornflowerblue", linewidth=2)
ax.plot(df["position"], df["alt_smooth"], label="ALT", color="firebrick", linewidth=2)

ax.axvline(197081044, color="gold", linewidth=1.5, linestyle="--")

# 4) Formatting

ax.set_xlim(197076044, 197086544)
ax.ticklabel_format(style="plain", axis="x")
ax.get_xaxis().get_major_formatter().set_useOffset(False)

ax.set_ylabel("RNA-seq signal")
ax.set_xlabel("Genomic position")

plt.tight_layout()
plt.savefig(out_png, dpi=300)
plt.show()
