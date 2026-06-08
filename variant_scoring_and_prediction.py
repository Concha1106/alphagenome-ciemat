#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 27 10:54:38 2026

@author: e6260
"""

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
from alphagenome_key import get_dna_model
import pandas as pd
from pathlib import Path

# -----------------------------
# 1) General configuration
# -----------------------------

OUTPUT_DIR = Path("results/2026-05-27-exp-all-output-types")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

dna_model = get_dna_model()

variant = genome.Variant(
    chromosome="chr3",
    position=197081044,
    reference_bases="TACTC",
    alternate_bases="T",
)

interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)


# -----------------------------
# 2) SCORE_VARIANT
# -----------------------------

selected_scorers = [
    variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITE_USAGE"],
    variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITES"],
    variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_JUNCTIONS"],
    variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"],
    variant_scorers.RECOMMENDED_VARIANT_SCORERS["ATAC"],
    variant_scorers.RECOMMENDED_VARIANT_SCORERS["DNASE"],
    variant_scorers.RECOMMENDED_VARIANT_SCORERS["CHIP_TF"],
    variant_scorers.RECOMMENDED_VARIANT_SCORERS["CHIP_HISTONE"],
    variant_scorers.RECOMMENDED_VARIANT_SCORERS["CONTACT_MAPS"],
    variant_scorers.RECOMMENDED_VARIANT_SCORERS["CAGE"],
    variant_scorers.RECOMMENDED_VARIANT_SCORERS["PROCAP"],
]

scores = dna_model.score_variant(
    interval=interval,
    variant=variant,   
    variant_scorers= selected_scorers,
)

scores_df = variant_scorers.tidy_scores(scores)
scores_df.head()


scores_df.to_csv(
    OUTPUT_DIR / "score-variant-all-output-types.tsv",
    sep="\t",
    index=False
)

print(scores_df.head())
print(scores_df.columns)

#dlg1_tibial = scores_df[
#    (scores_df["gene_name"] == "DLG1") &
#    (scores_df["ontology_curie"] == "UBERON:0007610")
#]

#dlg1_tibial.shape

# -----------------------------
# 3) BASIC SUMMARY TABLES
# -----------------------------

print("\nNumber of rows by output type:")
print(scores_df["output_type"].value_counts())

print("\nTop scores by quantile score:")
top_quantile = scores_df.sort_values(
    "quantile_score",
    ascending=False
).head(20)

print(top_quantile[
    [
        "output_type",
        "variant_scorer",
        "gene_name",
        "track_name",
        "ontology_curie",
        "biosample_name",
        "gtex_tissue",
        "raw_score",
        "quantile_score",
    ]
])

top_quantile.to_csv(
    OUTPUT_DIR / "score-variant-top-quantile-scores.tsv",
    sep="\t",
    index=False
)

dlg1_tibial = scores_df[
    (scores_df["gene_name"] == "DLG1")
    &
    (
        (scores_df["gtex_tissue"] == "Artery_Tibial")
        |
        (scores_df["ontology_curie"] == "UBERON:0007610")
    )
]

scores_df[
    scores_df["output_type"] == "SPLICE_SITES"
][
    [
        "gene_name",
        "ontology_curie",
        "gtex_tissue",
        "raw_score",
        "quantile_score"
    ]
]




requested_outputs = [
    dna_client.OutputType.RNA_SEQ,
    dna_client.OutputType.CAGE,
    dna_client.OutputType.PROCAP,
    dna_client.OutputType.SPLICE_SITES,
    dna_client.OutputType.SPLICE_SITE_USAGE,
    dna_client.OutputType.SPLICE_JUNCTIONS,
    dna_client.OutputType.DNASE,
    dna_client.OutputType.ATAC,
    dna_client.OutputType.CHIP_HISTONE,
    dna_client.OutputType.CHIP_TF,
    dna_client.OutputType.CONTACT_MAPS,
]