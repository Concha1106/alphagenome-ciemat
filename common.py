#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 14:41:38 2026

@author: e6260
"""

def classify_prediction_output(data):
    """
    Classify one AlphaGenome predict_variant output by its actual object structure.
    """

    if data is None:
        return "missing"

    if hasattr(data, "junctions"):
        return "junction_data"

    if hasattr(data, "values") and data.values.ndim == 3:
        return "contact_map_data"

    if hasattr(data, "values") and data.values.ndim == 2:
        return "track_data"

    return "unsupported"

