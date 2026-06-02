#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 13:57:54 2026

Script: Variant splicing scoring using AlphaGenome

Description:
This script evaluates the splicing impact of a genomic variant using AlphaGenome.
It computes splice site usage, splice site strength, and splice junction effects,
and integrates them into a combined splicing score.

@author: concha
"""

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
from alphagenome_key import get_dna_model
import pandas as pd

# 1) Define paper variable: 4pb deletion

variant = genome.Variant (
    chromosome="chr3",
    position= 197081044,
    reference_bases="TACTC",
    alternate_bases="T",
    )

# 2) Create model input interval

interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

# 3) Create AlphaGenome client
model = get_dna_model()

# 4) Define scores for splicing-related parameters
scores = model.score_variant(
    interval=interval,
    variant=variant,
    variant_scorers=[
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITE_USAGE"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITES"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_JUNCTIONS"],
    ],
)


# 5) Separate results
ad_usage = scores[0]
ad_sites = scores[1]
ad_junctions = scores[2]

# A) SPLICE_SITE_USAGE

df_usage = pd.DataFrame(
    ad_usage.X,
    index=ad_usage.obs["gene_name"],
    columns=ad_usage.var["gtex_tissue"],
)

df_usage_tibial = df_usage.loc[:, df_usage.columns.str.contains("Tibial", na=False)]
df_usage_tibial_clean = df_usage_tibial.T.groupby(level=0).mean().T


# B) SPLICE_SITES

df_sites = pd.DataFrame(
    ad_sites.X,
    index=ad_sites.obs["gene_name"],
    columns=ad_sites.var["name"] + "_" + ad_sites.var["strand"],
)


# C) SPLICE_JUNCTIONS
# ad_junctions.X contains junction scores per tissue track

ad_junctions.obs["length"] = (
    ad_junctions.obs["junction_End"] - ad_junctions.obs["junction_Start"]
)

df_junctions = pd.DataFrame(
    ad_junctions.X,
    index=ad_junctions.obs.index,
    columns=ad_junctions.var["gtex_tissue"],
)

df_junctions_tibial = df_junctions.loc[:, df_junctions.columns.str.contains("Tibial", na=False)].copy()
df_junctions_tibial["junction_id"] = df_junctions_tibial.index
df_junctions_tibial["junction_Start"] = ad_junctions.obs["junction_Start"].values
df_junctions_tibial["junction_End"] = ad_junctions.obs["junction_End"].values
df_junctions_tibial["length"] = ad_junctions.obs["length"].values
df_junctions_tibial["max_tibial_score"] = (
    df_junctions_tibial[["Artery_Tibial", "Nerve_Tibial"]].abs().max(axis=1)
)

# D) Merge splicing score

max_sites = float(df_sites.abs().max(axis=1).iloc[0])
max_usage = float(df_usage_tibial_clean.abs().max(axis=1).iloc[0])
max_junctions = float(
    df_junctions_tibial[["Artery_Tibial", "Nerve_Tibial"]].abs().to_numpy().max()
)

combined_score = max_sites + max_usage + (max_junctions / 5)

# E) Summary table

summary = pd.DataFrame({
    "gene_name": ["DLG1"],
    "usage_Artery_Tibial": [df_usage_tibial_clean.loc["DLG1", "Artery_Tibial"]],
    "usage_Nerve_Tibial": [df_usage_tibial_clean.loc["DLG1", "Nerve_Tibial"]],
    "sites_donor_minus": [df_sites.loc["DLG1", "donor_-"]],
    "sites_acceptor_minus": [df_sites.loc["DLG1", "acceptor_-"]],
})



# Longest junction (candidate for exon skipping)
longest_idx = df_junctions_tibial["length"].idxmax()
summary["candidate_junction_id"] = longest_idx
summary["candidate_junction_start"] = df_junctions_tibial.loc[longest_idx, "junction_Start"]
summary["candidate_junction_end"] = df_junctions_tibial.loc[longest_idx, "junction_End"]
summary["candidate_junction_length"] = df_junctions_tibial.loc[longest_idx, "length"]
summary["candidate_junction_Artery_Tibial"] = df_junctions_tibial.loc[longest_idx, "Artery_Tibial"]
summary["candidate_junction_Nerve_Tibial"] = df_junctions_tibial.loc[longest_idx, "Nerve_Tibial"]
summary["candidate_junction_max_tibial_score"] = df_junctions_tibial.loc[longest_idx, "max_tibial_score"]
summary["max_sites"] = max_sites
summary["max_usage"] = max_usage
summary["max_junctions"] = max_junctions
summary["combined_splicing_score"] = combined_score


# F) Save results

df_junctions_tibial = df_junctions_tibial.sort_values(
    by="max_tibial_score",
    ascending=False
)

outdir = "~/Desktop/alphagenome-ciemat/results/2026-04-10-exp02-exon-skipping"

df_usage_tibial_clean.to_csv(f"{outdir}/dlg1_usage_tibial.tsv", sep="\t")
df_sites.to_csv(f"{outdir}/dlg1_sites.tsv", sep="\t")
df_junctions_tibial.to_csv(f"{outdir}/dlg1_junctions_tibial.tsv", sep="\t", index=False)
summary.to_csv(f"{outdir}/dlg1_summary.tsv", sep="\t", index=False)

print("Results saved in:")
print(f"- {outdir}/dlg1_usage_tibial.tsv")
print(f"- {outdir}/dlg1_sites.tsv")
print(f"- {outdir}/dlg1_junctions_tibial.tsv")
print("Combined splicing score =", combined_score)
print(f"- {outdir}/dlg1_summary.tsv")
print("Candidate junction ID =", summary.loc[0, "candidate_junction_id"])
print(
    "Candidate junction =",
    f"{summary.loc[0, 'candidate_junction_start']} -> "
    f"{summary.loc[0, 'candidate_junction_end']}"
)
print("Candidate junction length =", summary.loc[0, "candidate_junction_length"])
print("Candidate junction max tibial score =", summary.loc[0, "candidate_junction_max_tibial_score"])






