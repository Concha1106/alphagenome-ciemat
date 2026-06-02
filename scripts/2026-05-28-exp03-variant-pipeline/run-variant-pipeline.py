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
import pandas as pd

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
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
    parser.add_argument("--ontology-curie", required=True, help="Ontology CURIE used for predict_variant, e.g. UBERON:0007610")
    
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
# 4) Alphagenome model call and Output configuration
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

def get_ontology_terms(ontology_curie):
    """
    Define ontology terms for predict_variant.

    In the current MVP, one ontology CURIE is required to avoid requesting
    predictions for all available ontologies.
    """

    return [ontology_curie]

# -----------------------------------------------
# 5) Score variant analysis
# -----------------------------------------------

def run_score_variant(dna_model, interval, variant):
    """
    Run AlphaGenome score_variant() using recommended variant scorers.

    If no variant_scorers are provided, AlphaGenome uses the recommended
    scorers for the selected organism.
    """

    scores = dna_model.score_variant(
        interval=interval,
        variant=variant
    )

    return scores


def tidy_score_results(scores):
    """
    Convert AlphaGenome score_variant output into a tidy pandas DataFrame.
    """

    scores_df = variant_scorers.tidy_scores(scores)

    return scores_df


def export_score_tables(scores_df, output_dir):
    """
    Export complete score table.

    No prioritization is applied at this stage.
    The goal is to preserve all information generated by score_variant().
    """

    all_scores_path = output_dir / "score_variant_all.tsv"

    scores_df.to_csv(
        all_scores_path,
        sep="\t",
        index=False
    )

    return all_scores_path

# -----------------------------------------------
# 6) Predict variant analysis
# -----------------------------------------------

def run_predict_variant(dna_model, interval, variant, requested_outputs, ontology_terms):
    """
    Run AlphaGenome predict_variant() for REF and ALT alleles.

    In the current MVP, ontology_terms must be provided to avoid requesting
    all available ontology-specific predictions.
    """

    prediction = dna_model.predict_variant(
        interval=interval,
        variant=variant,
        requested_outputs=requested_outputs,
        ontology_terms=ontology_terms
    )

    return prediction

def export_one_trackdata_output(prediction, output_name, output_dir):
    """
    Export REF, ALT and delta values for one TrackData output.

    This first version is intended to validate the export logic using one
    output, for example rna_seq.
    """

    ref_data = getattr(prediction.reference, output_name)
    alt_data = getattr(prediction.alternate, output_name)

    values_ref = ref_data.values
    values_alt = alt_data.values

    metadata = ref_data.metadata.copy()

    rows = []

    for track_index in range(values_ref.shape[1]):
        track_metadata = metadata.iloc[track_index].to_dict()

        for position_index in range(values_ref.shape[0]):
            ref_value = values_ref[position_index, track_index]
            alt_value = values_alt[position_index, track_index]

            genomic_position_0based = ref_data.interval.start + position_index
            genomic_position_1based = genomic_position_0based + 1

            row = {
                "output_name": output_name,
                "position_index": position_index,
                "genomic_position_0based": genomic_position_0based,
                "genomic_position_1based": genomic_position_1based,
                "chromosome": ref_data.interval.chromosome,
                "track_index": track_index,
                "ref_value": ref_value,
                "alt_value": alt_value,
                "delta_alt_ref": alt_value - ref_value,
                }
            row.update(track_metadata)

            rows.append(row)

    output_df = pd.DataFrame(rows)

    output_path = output_dir / f"predict_variant_{output_name}.tsv"

    output_df.to_csv(
        output_path,
        sep="\t",
        index=False
    )

    return output_path


def export_all_trackdata_outputs(prediction, output_dir):
    """
    Export REF, ALT and delta values for all TrackData-like outputs.

    Outputs with zero tracks are skipped.
    """

    trackdata_outputs = [
        "atac",
        "cage",
        "chip_histone",
        "chip_tf",
        "dnase",
        "procap",
        "rna_seq",
        "splice_sites",
        "splice_site_usage",
    ]

    exported_paths = []

    for output_name in trackdata_outputs:
        ref_data = getattr(prediction.reference, output_name)

        if ref_data.values.shape[1] == 0:
            print(f"Skipping {output_name}: no tracks available.")
            continue

        output_path = export_one_trackdata_output(
            prediction,
            output_name,
            output_dir
        )

        exported_paths.append(output_path)

    return exported_paths

def export_splice_junctions(prediction, output_dir):
    """
    Export REF, ALT and delta values for splice_junctions.

    In JunctionData, rows correspond to predicted splice junctions and columns
    correspond to tracks.
    """

    ref_data = prediction.reference.splice_junctions
    alt_data = prediction.alternate.splice_junctions

    values_ref = ref_data.values
    values_alt = alt_data.values

    metadata = ref_data.metadata.copy()
    junctions = ref_data.junctions.copy()
    
    

    if values_ref.shape[1] == 0:
        print("Skipping splice_junctions: no tracks available.")
        return None

    rows = []

    for track_index in range(values_ref.shape[1]):
        track_metadata = metadata.iloc[track_index].to_dict()

        for junction_index in range(values_ref.shape[0]):
            ref_value = values_ref[junction_index, track_index]
            alt_value = values_alt[junction_index, track_index]

            row = {
                "output_name": "splice_junctions",
                "junction_index": junction_index,
                "track_index": track_index,
                "ref_value": ref_value,
                "alt_value": alt_value,
                "delta_alt_ref": alt_value - ref_value,
            }

            row.update(track_metadata)

            junction = junctions[junction_index]

            row.update({
                "junction_chromosome": junction.chromosome,
                "junction_start_0based": junction.start,
                "junction_end_0based": junction.end,
                "junction_start_1based": junction.start + 1,
                "junction_end_1based": junction.end,
                "junction_strand": junction.strand,
                "junction_name": junction.name,
                })

            rows.append(row)

    output_df = pd.DataFrame(rows)

    output_path = output_dir / "predict_variant_splice_junctions.tsv"

    output_df.to_csv(
        output_path,
        sep="\t",
        index=False
    )

    return output_path

def export_contact_maps(prediction, output_dir):
    """
    Export REF, ALT and delta values for contact_maps.

    Contact maps are represented as matrices:
    bin_1 x bin_2 x track.
    """

    ref_data = prediction.reference.contact_maps
    alt_data = prediction.alternate.contact_maps

    values_ref = ref_data.values
    values_alt = alt_data.values

    metadata = ref_data.metadata.copy()

    if values_ref.shape[2] == 0:
        print("Skipping contact_maps: no tracks available.")
        return None

    rows = []

    for track_index in range(values_ref.shape[2]):
        track_metadata = metadata.iloc[track_index].to_dict()

        for bin_1 in range(values_ref.shape[0]):
            for bin_2 in range(values_ref.shape[1]):
                ref_value = values_ref[bin_1, bin_2, track_index]
                alt_value = values_alt[bin_1, bin_2, track_index]

                row = {
                    "output_name": "contact_maps",
                    "bin_1_index": bin_1,
                    "bin_2_index": bin_2,
                    "track_index": track_index,
                    "ref_value": ref_value,
                    "alt_value": alt_value,
                    "delta_alt_ref": alt_value - ref_value,
                }

                row.update(track_metadata)
                rows.append(row)

    output_df = pd.DataFrame(rows)

    output_path = output_dir / "predict_variant_contact_maps.tsv"

    output_df.to_csv(
        output_path,
        sep="\t",
        index=False
    )

    return output_path

def export_global_delta_tops(output_dir, top_n):
    """
    Export global top positive and negative delta tables from predict_variant outputs.

    This function reads predict_variant TSV files already generated by the pipeline,
    combines them, and exports global top increased/decreased predictions.
    """

    prediction_files = list(output_dir.glob("predict_variant_*.tsv"))

    prediction_files = [
        path for path in prediction_files
        if not path.name.endswith("_top_positive.tsv")
        and not path.name.endswith("_top_negative.tsv")
    ]

    dataframes = []

    for path in prediction_files:
        df = pd.read_csv(path, sep="\t")

        if "delta_alt_ref" not in df.columns:
            continue

        dataframes.append(df)

    if not dataframes:
        print("No predict_variant tables with delta_alt_ref found.")
        return None, None

    combined_df = pd.concat(dataframes, ignore_index=True, sort=False)

    top_positive = (
    combined_df
    .sort_values(by="delta_alt_ref", ascending=False)
    .groupby("output_name", group_keys=False)
    .head(top_n)
    )

    top_negative = (
    combined_df
    .sort_values(by="delta_alt_ref", ascending=True)
    .groupby("output_name", group_keys=False)
    .head(top_n)
    )
    
    top_positive_path = output_dir / "predict_variant_top_positive.tsv"
    top_negative_path = output_dir / "predict_variant_top_negative.tsv"

    top_positive.to_csv(top_positive_path, sep="\t", index=False)
    top_negative.to_csv(top_negative_path, sep="\t", index=False)

    return top_positive_path, top_negative_path
    

# -----------------------------------------------
# 7) Run log generation
# -----------------------------------------------

def write_runlog(
    args,
    output_dir,
    variant,
    interval,
    requested_outputs,
    ontology_terms,
    score_table_path,
    trackdata_paths,
    splice_junctions_path,
    contact_maps_path,
    top_positive_path,
    top_negative_path,
):
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
        runlog.write(f"predict_variant_ontology_terms: {ontology_terms}\n")
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
        
        runlog.write("\nGenerated files\n")
        runlog.write("-" * 20 + "\n")
        runlog.write(f"{score_table_path.name}\n")

        for path in trackdata_paths:
            runlog.write(f"{path.name}\n")

        if splice_junctions_path is not None:
            runlog.write(f"{splice_junctions_path.name}\n")

        if contact_maps_path is not None:
            runlog.write(f"{contact_maps_path.name}\n")
        else:
            runlog.write("predict_variant_contact_maps.tsv: not generated, no tracks available\n")

        if top_positive_path is not None:
            runlog.write(f"{top_positive_path.name}\n")

        if top_negative_path is not None:
            runlog.write(f"{top_negative_path.name}\n")
        
# -----------------------------------------------
# 8) Main workflow
# -----------------------------------------------

def main():

    args = parse_arguments()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    variant = build_variant(args)
    interval = build_interval(variant, args.interval_size)
    
    dna_model = create_dna_model()
    requested_outputs = get_requested_outputs(args.output_types)
    
    scores = run_score_variant(dna_model, interval, variant)
    scores_df = tidy_score_results(scores)
    score_table_path = export_score_tables(scores_df, output_dir)
    
    ontology_terms = get_ontology_terms(args.ontology_curie)
    prediction = run_predict_variant(dna_model, interval, variant, requested_outputs, ontology_terms)
    trackdata_paths = export_all_trackdata_outputs(prediction, output_dir)
    splice_junctions_path = export_splice_junctions(prediction, output_dir)
    contact_maps_path = export_contact_maps(prediction, output_dir)
    top_positive_path, top_negative_path = export_global_delta_tops(output_dir, args.top_n)
    
    write_runlog(
    args,
    output_dir,
    variant,
    interval,
    requested_outputs,
    ontology_terms,
    score_table_path,
    trackdata_paths,
    splice_junctions_path,
    contact_maps_path,
    top_positive_path,
    top_negative_path,
)
    
    print("AlphaGenome DNA model client created successfully.")
    print("score_variant() completed successfully.")
    print(f"Complete score table saved to: {score_table_path}")
    print("predict_variant() completed successfully.")
    print(f"TrackData prediction tables saved: {len(trackdata_paths)}")
    print(f"Splice junctions prediction table saved to: {splice_junctions_path}")
    print(f"Contact maps prediction table saved to: {contact_maps_path}")
    print(f"Top positive delta table saved to: {top_positive_path}")
    print(f"Top negative delta table saved to: {top_negative_path}")
    print("Runlog created successfully.")
        



if __name__ == "__main__":
    main()