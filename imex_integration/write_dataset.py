"""
Module to rewrite specific keywords from IMEX data set after collecting upscaling informations of permeability, porosity and primal mesh
"""
import numpy as np

def __init__(self, original_dataset, new_dataset):

    porosity_token = '*POR'
    permeability_token = '*PERM'
    well_token = '*WELL'

    with open(original_dataset, 'r') as dataset:
        original_dataset = dataset.readlines()

def write_effective_porosity(self):
    pass

def write_effective_permeability(self):
    pass

def well_settings(self):
    pass
