#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 16:10:52 2026
Script title: General variant evaluation pipeline - initial version
Description:
    Initial script for experiment exp03.
    This version reads a genomic variant from command-line arguments,
    creates the output directory and writes a basic runlog.

@author: concha
"""
import argparse
from datetime import datetime
from pathlib import Path


# ---------------------------------
# 1) Argument parsing
# ---------------------------------

def parse_arguments():
    """
    Define and read command-line arguments for the variant evaluation pipeline.
    """

    parser = argparse.ArgumentParser(
        description="General variant evaluation pipeline with AlphaGenome."
    )
    # Future development:
    # - quantile-based prioritization of score_variant() results.
    # - predict_mode = all | prioritized.
    # These features are planned but not implemented in the current version.


    # Mandatory arguments: define the genomic variant.
    parser.add_argument("--chrom", required=True, help="Chromosome, e.g. chr3")
    parser.add_argument("--pos", required=True, type=int, help="Genomic position, 1-based")
    parser.add_argument("--ref", required=True, help="Reference allele")
    parser.add_argument("--alt", required=True, help="Alternative allele")

    # General execution options.
    parser.add_argument(
        "--interval-size",
        default="1MB",
        choices=["1MB", "500KB", "100KB", "16KB"],
        help="Sequence interval size used for AlphaGenome analysis"
    )

    parser.add_argument(
        "--output-dir",
        default="results/2026-05-28-exp03-variant-pipeline",
        help="Directory where output files will be saved"
    )

    parser.add_argument(
        "--top-n",
        default=10,
        type=int,
        help="Number of top results to export in summary tables"
    )


    # Optional metadata filters: These filters generate extra views only.
    # They must not removeinformation from complete output tables

    parser.add_argument("--ontology-curie", default=None, help="Optional ontology CURIE filter")
    parser.add_argument("--gtex-tissue", default=None, help="Optional GTEx tissue filter")
    parser.add_argument("--biosample-name", default=None, help="Optional biosample name filter")
    parser.add_argument("--gene-name", default=None, help="Optional gene name filter")
    parser.add_argument("--gene-id", default=None, help="Optional Ensembl gene ID filter")

    parser.add_argument(
        "--output-types",
        nargs="+",
        default=None,
        help="Optional list of AlphaGenome output types to analyse with predict_variant"
    )

    return parser.parse_args()


# ---------------------------------
# 2) Run log generation
# ---------------------------------

def write_runlog(args, output_dir):
    """
    Write a basic runlog for the current execution.

    The runlog records:
    - execution date;
    - input variant;
    - execution parameters;
    - optional filters;
    - output directory.
    """

    runlog_path = output_dir / "runlog.txt"

    with open(runlog_path, "w", encoding="utf-8") as runlog:
        runlog.write("Experiment exp03 — General Variant Evaluation Pipeline\n")
        runlog.write("=" * 60 + "\n\n")

        runlog.write(f"Execution date: {datetime.now()}\n\n")

        runlog.write("Input variant\n")
        runlog.write("-" * 20 + "\n")
        runlog.write(f"chromosome: {args.chrom}\n")
        runlog.write(f"position_1based: {args.pos}\n")
        runlog.write(f"reference: {args.ref}\n")
        runlog.write(f"alternate: {args.alt}\n\n")

        runlog.write("Execution parameters\n")
        runlog.write("-" * 20 + "\n")
        runlog.write(f"interval_size: {args.interval_size}\n")
        runlog.write(f"top_n: {args.top_n}\n")
        runlog.write(f"output_types: {args.output_types}\n\n")

        runlog.write("Optional secondary-view filters\n")
        runlog.write("-" * 20 + "\n")
        runlog.write(f"ontology_curie: {args.ontology_curie}\n")
        runlog.write(f"gtex_tissue: {args.gtex_tissue}\n")
        runlog.write(f"biosample_name: {args.biosample_name}\n")
        runlog.write(f"gene_name: {args.gene_name}\n")
        runlog.write(f"gene_id: {args.gene_id}\n\n")

        runlog.write("Output directory\n")
        runlog.write("-" * 20 + "\n")
        runlog.write(f"{output_dir}\n")



# ---------------------------------
# 3) Main workflow
# ---------------------------------

def main():

    args = parse_arguments()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    write_runlog(args, output_dir)

    print("Runlog created successfully.")


if __name__ == "__main__":
    main()