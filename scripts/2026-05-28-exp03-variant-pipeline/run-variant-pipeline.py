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
# -----------------------------------------------
# 1) Imports
# -----------------------------------------------

import argparse
from datetime import datetime
from pathlib import Path

from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome_key import get_dna_model


# -----------------------------------------------
# 2) Argument parsing
# -----------------------------------------------


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
    # They must not remove information from complete output tables

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



# -----------------------------------------------
# 3) Variant validation and interval definition 
# -----------------------------------------------

def normalize_chromosome(chrom):
    """
    Normalize chromosome format.

    Example:
    - "3" becomes "chr3"
    - "chr3" remains "chr3"
    """

    chrom = chrom.strip()

    if not chrom.startswith("chr"):
        chrom = "chr" + chrom

    return chrom


def normalize_allele(allele):
    """
    Normalize REF or ALT allele.

    The allele is converted to uppercase and surrounding spaces are removed.
    """

    allele = allele.strip().upper()

    return allele

def validate_position(pos):
    """
    Validate genomic position.

    The input position must be a positive integer.
    AlphaGenome Variant.position uses 1-based coordinates.
    """

    if pos < 1:
        raise ValueError("Genomic position must be a positive 1-based integer.")

    return pos

def build_variant(args):
    """
    Build an AlphaGenome Variant object from command-line arguments.

    The Variant object is used as the main internal representation of the
    genomic variant in the pipeline.
    """

    chrom = normalize_chromosome(args.chrom)
    pos = validate_position(args.pos)
    ref = normalize_allele(args.ref)
    alt = normalize_allele(args.alt)

    variant = genome.Variant(
        chromosome=chrom,
        position=pos,
        reference_bases=ref,
        alternate_bases=alt
    )

    return variant

def get_sequence_length(interval_size):
    """
    Convert interval size argument into an AlphaGenome sequence length.
    """

    sequence_length_map = {
        "1MB": dna_client.SEQUENCE_LENGTH_1MB,
        "500KB": dna_client.SEQUENCE_LENGTH_500KB,
        "100KB": dna_client.SEQUENCE_LENGTH_100KB,
        "16KB": dna_client.SEQUENCE_LENGTH_16KB,
    }

    return sequence_length_map[interval_size]

def build_interval(variant, interval_size):
    """
    Build AlphaGenome-compatible interval around the variant.
    """

    sequence_length = get_sequence_length(interval_size)

    interval = variant.reference_interval.resize(
        sequence_length
    )

    return interval

# -----------------------------------------------
# 3) Alphagenome model call and Output configuration
# -----------------------------------------------

def create_dna_model():
    """
    Create AlphaGenome DNA model client using the project helper function.

    The API key is read from the ALPHAGENOME_API_KEY environment variable
    inside alphagenome_key.py.
    """

    dna_model = get_dna_model()

    return dna_model

def get_requested_outputs(output_types):
    """
    Define AlphaGenome output types requested for predict_variant.

    If output_types is None, all AlphaGenome output types are used.
    """

    if output_types is None:
        return list(dna_client.OutputType)

    requested_outputs = []

    for output_type in output_types:
        requested_outputs.append(dna_client.OutputType[output_type])

    return requested_outputs


# -----------------------------------------------
# X) Run log generation
# -----------------------------------------------

def write_runlog(args, output_dir, variant, interval, requested_outputs):
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
        
        runlog.write("\nVariant classification\n")
        runlog.write("-" * 20 + "\n")
        runlog.write(f"is_snv: {variant.is_snv}\n")
        runlog.write(f"is_insertion: {variant.is_insertion}\n")
        runlog.write(f"is_deletion: {variant.is_deletion}\n")
        runlog.write(f"is_indel: {variant.is_indel}\n")
        runlog.write(f"is_structural: {variant.is_structural}\n")
        runlog.write(f"is_frameshift: {variant.is_frameshift}\n")
        runlog.write(f"variant_start_0based: {variant.start}\n")
        runlog.write(f"variant_end_0based: {variant.end}\n")
        
        runlog.write("\nInterval definition\n")
        runlog.write("-" * 20 + "\n")
        runlog.write(f"reference_interval: {variant.reference_interval}\n")
        runlog.write(f"input_interval: {interval}\n")
        runlog.write(f"interval_chromosome: {interval.chromosome}\n")
        runlog.write(f"interval_start_0based: {interval.start}\n")
        runlog.write(f"interval_end_0based: {interval.end}\n")
        runlog.write(f"interval_width: {interval.width}\n")

        runlog.write("\nAlphaGenome model configuration\n")
        runlog.write("-" * 20 + "\n")
        runlog.write("model_client: created with get_dna_model()\n")
        runlog.write("api_key_source: ALPHAGENOME_API_KEY environment variable\n")
        runlog.write(f"number_requested_outputs: {len(requested_outputs)}\n")
        runlog.write(
            "requested_outputs: "
            + ", ".join([output.name for output in requested_outputs])
            + "\n"
        )

# -----------------------------------------------
# X) Main workflow
# -----------------------------------------------

def main():

    args = parse_arguments()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    variant = build_variant(args)
    interval = build_interval(variant, args.interval_size)
    
    dna_model = create_dna_model()
    requested_outputs = get_requested_outputs(args.output_types)

    
    write_runlog(args, output_dir, variant, interval, requested_outputs)
    
    print("AlphaGenome DNA model client created successfully.")
    print("Runlog created successfully.")
    
    


if __name__ == "__main__":
    main()