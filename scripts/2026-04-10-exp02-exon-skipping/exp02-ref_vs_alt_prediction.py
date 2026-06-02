#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script: exp02-ref_vs_alt_prediction.py

Description:
Generate REF vs ALT AlphaGenome predictions for the DLG1 splicing variant
reported in the paper (chr3:197081044 TACTC>T), focusing on outputs needed
to reproduce an exon-skipping-like visualization:

- splice junctions
- splice site usage
- splice sites
- RNA-seq predicted coverage

The output tables are intended for downstream plotting.

@author: concha
"""

from pathlib import Path
import pandas as pd
from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.models.dna_output import OutputType

from alphagenome_key import get_dna_model



# 1) Define variant and prediction interval

variant = genome.Variant(
    chromosome="chr3",
    position=197081044,
    reference_bases="TACTC",
    alternate_bases="T",
)

interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)


# 2) Load AlphaGenome model

model = get_dna_model()


# 3) Run REF vs ALT prediction

variant_output = model.predict_variant(
    interval=interval,
    variant=variant,
    requested_outputs=[
        OutputType.SPLICE_SITES,
        OutputType.SPLICE_SITE_USAGE,
        OutputType.SPLICE_JUNCTIONS,
        OutputType.RNA_SEQ,
    ],
    ontology_terms=["UBERON:0007610"],  # GTEx Artery_Tibial
)


# 4) Define output directory and plotting region

outdir = Path("~/Desktop/test-zip-recovery/alphagenome-ciemat/results/2026-04-10-exp02-exon-skipping").expanduser()
outdir.mkdir(parents=True, exist_ok=True)

plot_region = genome.Interval(
    chromosome="chr3",
    start=197076044,
    end  =197086544,
)

print(variant_output.reference.splice_junctions.metadata)

# 5) Helper function for TrackData outputs

def track_to_dataframe(ref_track, alt_track, region, output_type):
    # Convert REF and ALT TrackData objects into a tidy dataframe.

    ref_region = ref_track.slice_by_interval(region)
    alt_region = alt_track.slice_by_interval(region)

    positions = list(range(region.start, region.end))
    dfs = []

    for i, track_name in enumerate(ref_region.names):
        df = pd.DataFrame({
            "chromosome": region.chromosome,
            "position": positions,
            "output_type": output_type,
            "track": str(track_name),
            "track_index": i,
            "ref_value": ref_region.values[:, i],
            "alt_value": alt_region.values[:, i],
        })

        df["delta_alt_ref"] = df["alt_value"] - df["ref_value"]
        dfs.append(df)

    return pd.concat(dfs, ignore_index=True)


# 6) Splice junctions: REF vs ALT

ref_junc = variant_output.reference.splice_junctions.filter_to_negative_strand()
alt_junc = variant_output.alternate.splice_junctions.filter_to_negative_strand()

df_junc = pd.DataFrame({
    "chromosome": [j.chromosome for j in ref_junc.junctions],
    "start": [j.start for j in ref_junc.junctions],
    "end": [j.end for j in ref_junc.junctions],
    "strand": [j.strand for j in ref_junc.junctions],
    "ref_value": ref_junc.values[:, 0],
    "alt_value": alt_junc.values[:, 0],
})

df_junc["delta_alt_ref"] = df_junc["alt_value"] - df_junc["ref_value"]

df_junc.to_csv(
    outdir / "dlg1_splice_junctions_ref_vs_alt.tsv",
    sep="\t",
    index=False,
)

df_junc_top_up = df_junc.sort_values("delta_alt_ref", ascending=False).head(10)
df_junc_top_down = df_junc.sort_values("delta_alt_ref", ascending=True).head(10)

df_junc_top_up.to_csv(
    outdir / "dlg1_splice_junctions_top10_increased_alt_vs_ref.tsv",
    sep="\t",
    index=False,
)

df_junc_top_down.to_csv(
    outdir / "dlg1_splice_junctions_top10_decreased_alt_vs_ref.tsv",
    sep="\t",
    index=False,
)

# 7) Splice site usage: REF vs ALT

ref_usage = variant_output.reference.splice_site_usage.filter_to_negative_strand()
alt_usage = variant_output.alternate.splice_site_usage.filter_to_negative_strand()

df_usage = track_to_dataframe(
    ref_usage,
    alt_usage,
    plot_region,
    "splice_site_usage",
)

df_usage.to_csv(
    outdir / "dlg1_splice_site_usage_ref_vs_alt.tsv",
    sep="\t",
    index=False,
)

df_usage_top_up = df_usage.sort_values("delta_alt_ref", ascending=False).head(10)
df_usage_top_down = df_usage.sort_values("delta_alt_ref", ascending=True).head(10)

df_usage_top_up.to_csv(
    outdir / "dlg1_splice_site_usage_top10_increased_alt_vs_ref.tsv",
    sep="\t",
    index=False,
)

df_usage_top_down.to_csv(
    outdir / "dlg1_splice_site_usage_top10_decreased_alt_vs_ref.tsv",
    sep="\t",
    index=False,
)
# 8) Splice sites: REF vs ALT

ref_sites = variant_output.reference.splice_sites.filter_to_negative_strand()
alt_sites = variant_output.alternate.splice_sites.filter_to_negative_strand()

df_sites = track_to_dataframe(
    ref_sites,
    alt_sites,
    plot_region,
    "splice_sites",
)

df_sites.to_csv(
    outdir / "dlg1_splice_sites_ref_vs_alt.tsv",
    sep="\t",
    index=False,
)

df_sites_top_up = df_sites.sort_values("delta_alt_ref", ascending=False).head(10)
df_sites_top_down = df_sites.sort_values("delta_alt_ref", ascending=True).head(10)

df_sites_top_up.to_csv(
    outdir / "dlg1_splice_site_top10_increased_alt_vs_ref.tsv",
    sep="\t",
    index=False,
)

df_sites_top_down.to_csv(
    outdir / "dlg1_splice_site_top10_decreased_alt_vs_ref.tsv",
    sep="\t",
    index=False,
)

# 9) RNA-seq predicted coverage: REF vs ALT

df_rna = track_to_dataframe(
    variant_output.reference.rna_seq,
    variant_output.alternate.rna_seq,
    plot_region,
    "rna_seq",
)

df_rna.to_csv(
    outdir / "dlg1_rna_seq_ref_vs_alt.tsv",
    sep="\t",
    index=False,
)


# 10) Report outputs

print("Saved output tables in:", outdir)
print("- dlg1_splice_junctions_ref_vs_alt.tsv")
print("- dlg1_splice_site_usage_ref_vs_alt.tsv")
print("- dlg1_splice_sites_ref_vs_alt.tsv")
print("- dlg1_rna_seq_ref_vs_alt.tsv")
print("- dlg1_splice_junctions_top10_increased_alt_vs_ref.tsv")
print("- dlg1_splice_junctions_top10_decreased_alt_vs_ref.tsv")
print("- dlg1_splice_site_usage_top10_increased_alt_vs_ref.tsv")
print("- dlg1_splice_site_usage_top10_decreased_alt_vs_ref.tsv")
print("- dlg1_splice_site_top10_increased_alt_vs_ref.tsv")
print("- dlg1_splice_site_top10_decreased_alt_vs_ref.tsv")