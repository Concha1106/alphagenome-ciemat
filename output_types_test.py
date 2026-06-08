#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 14:38:30 2026

@author: e6260
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
        OutputType.ATAC,
        OutputType.DNASE,
        OutputType.CHIP_TF,
        OutputType.CHIP_HISTONE,
        OutputType.CAGE,
        OutputType.PROCAP,
        #OutputType.CONTACT_MAPS,
    ],
    ontology_terms=["UBERON:0007610"],  # GTEx Artery_Tibial
)