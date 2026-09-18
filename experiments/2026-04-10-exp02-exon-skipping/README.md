# Experiment exp02 — Exon Skipping (DLG1)

## Objective

Validate the analysis workflow by reproducing the DLG1 exon-skipping case described by Avsec et al. (2026), using both the effect metrics returned by `score_variant()` and the direct REF/ALT predictions generated with `predict_variant()`.

## Input data

- Variant: `chr3:197081044 TACTC>T`
- Gene: `DLG1`
- Reference genome: GRCh38
- Ontology: `UBERON:0007610` (Artery Tibial)
- MANE Select transcript: `ENST00000667157`
- Gene annotation: Ensembl GRCh38, release 115

The reference case corresponds to the exon-skipping event shown in Figure 3b of Avsec et al. (2026).

## Analysis

The experiment combines two AlphaGenome approaches:

1. `score_variant()`, used to obtain summary metrics for the predicted effect of the variant on splicing.
2. `predict_variant()`, used to compare REF and ALT predictions directly for:
   - RNA-seq
   - splice junctions
   - splice sites
   - splice site usage

Splicing-related signals were analyzed on the negative strand, corresponding to `DLG1`.

The MANE Select transcript annotation was used to interpret the predicted splice junctions in the context of the exon structure of `DLG1`.

## Scripts

- `scripts/2026-04-10-exp02-exon-skipping/exp02-score_exon_skipping.py`
- `scripts/2026-04-10-exp02-exon-skipping/exp02-ref_vs_alt_prediction.py`
- `scripts/2026-04-10-exp02-exon-skipping/plot-fig3b-junctions.py`
- `scripts/2026-04-10-exp02-exon-skipping/plot-fig3b-rna-seq.py`
- `scripts/2026-04-10-exp02-exon-skipping/plot-fig3b-splice-sites-and-usage.py`

## Main results

The `score_variant()` analysis indicated a relevant effect on splicing, with a combined score of approximately 3.61.

The REF/ALT comparison generated with `predict_variant()` reproduced the expected exon-skipping pattern:

- loss of the junctions connecting exons 18–17 and 17–16;
- appearance of an ALT junction connecting exons 18–16 directly;
- reduced RNA-seq signal over exon 17;
- concordant changes in splice sites and splice site usage.

Together, the predictions reproduce the exon 17 skipping event described by Avsec et al. (2026).

## Visualization

The main outputs are represented using dedicated scripts for RNA-seq, splice junctions, splice sites and splice site usage.

The visualizations use a shared genomic region and the MANE Select transcript annotation to support REF/ALT comparison and structural interpretation of the event.

## Associated files

**Intermediate data:** `interm-data/2026-04-10-exp02-exon-skipping/`

**Results:** `results/2026-04-10-exp02-exon-skipping/`

The `runlog.txt` file records the main execution parameters.

## Status

Experiment completed and used as a validation case for the analysis workflow.
