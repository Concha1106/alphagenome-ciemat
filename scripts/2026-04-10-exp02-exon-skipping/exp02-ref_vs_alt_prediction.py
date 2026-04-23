#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 17:30:39 2026

Script: TBD

Description: Generate REF vs ALT AlphaGenome predictions for the DLG1 splicing variant
reported in the paper (chr3:197081044 TACTC>T), focusing on splicing-related
outputs and RNA-seq in tibial artery-related tissue context

@author: concha
"""

from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.models.dna_output import OutputType
from alphagenome_key import get_dna_model

# 1) Define paper variant: 4-bp deletion in DLG1
variant = genome.Variant(
    chromosome="chr3",
    position=197081044,
    reference_bases="TACTC",
    alternate_bases="T",
)

# 2) Define model input interval (1 Mb centered on variant reference interval)
interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

# 3) Load AlphaGenome client
model = get_dna_model()

# 4) Selection of outputs relevant to exon-skipping interpretation


variant_output = model.predict_variant(
    interval=interval,
    variant=variant,
    requested_outputs=[
        OutputType.SPLICE_SITES,
        OutputType.SPLICE_SITE_USAGE,
        OutputType.SPLICE_JUNCTIONS,
        OutputType.RNA_SEQ,
    ],
    ontology_terms=["UBERON:0007610"],
)

