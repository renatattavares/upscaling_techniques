"""
Module to rewrite specific keywords from IMEX data set after collecting upscaling informations of permeability, porosity and primal mesh
"""
import numpy as np

class DatasetWriter:

    def __init__(self, original_dataset):

        porosity_token = '*POR'
        permeability_token = '*PERM'
        well_token = '*WELL'

        with open(original_dataset, 'r') as old:
            self.original_dataset = old.readlines()

        self.identify_heading()

        with open('coarse_model.dat', 'w') as new:
            new.writelines(self.heading)

    def identify_heading(self):

        heading_token = 'Reservoir Description Section'

        for line in self.original_dataset:
            if heading_token in line:
                index = self.original_dataset.index(line)

        self.heading = self.original_dataset[:index]
        separator = self.original_dataset[index-1]
        self.heading.append(heading_token + '\n')
        self.heading.append(separator)

    def effective_porosity(self):
        pass

    def effective_permeability(self):
        pass
