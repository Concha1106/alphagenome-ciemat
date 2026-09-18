# Experiment exp01 — Installation Test

## Objective

Verify the AlphaGenome installation and basic functionality, and establish the initial structure for reproducible project organization.

## Input data

- Ensembl GTF
  - Source: Ensembl
  - File: `Homo_sapiens.GRCh38.115.abinitio.gtf.gz`

- Simulated GTF
  - Source: test file generated for this project
  - File: `mini_clean2.gtf`

The input data are located in `raw-data/2026-04-06-exp01-installation-test/`.

## Scripts

- `scripts/2026-04-06-exp01-installation-test/process_gtf_test.py`
- `scripts/2026-04-06-exp01-installation-test/process_gtf.py`

## Results

- Successful execution of the initial processing tests.
- Generation of `mini_output.feather`.

## Traceability

External data used in the project are recorded in `docs/manifest-data.tsv`.

## Status

Experiment completed.
