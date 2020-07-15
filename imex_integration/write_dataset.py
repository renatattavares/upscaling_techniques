"""
Module to rewrite specific keywords from IMEX data set after collecting upscaling informations of permeability, porosity and primal mesh
"""
import numpy as np

class DatasetWriter:

    def __init__(self, original_dataset, coarsening, effective_porosity, effective_permeability):

        self.separator = '********************************************************************************'
        self.directions = np.array(['I', 'J', 'K'])

        with open(original_dataset, 'r') as old:
            self.original_dataset = old.readlines()

        with open('coarse_model.dat', 'w') as new:
            pass

        self.write_heading()
        self.write_mesh_settings(coarsening)

    def write_heading(self):

        heading_token = '**         Reservoir Description Section                                      **'

        for line in self.original_dataset:
            if heading_token in line:
                index = self.original_dataset.index(line)

        heading = self.original_dataset[:index]
        heading.append(heading_token + '\n')
        heading.append(self.separator + '\n')

        self.write_content(heading)

    def write_mesh_settings(self, coarsening):

        grid_token = '*GRID *CART'
        length_token = '*D *CON'
        mesh_settings = grid_token

        for dir in coarsening:
            mesh_settings = mesh_settings + ' ' + str(dir)

        mesh_settings = mesh_settings + '\n*KDIR *KDOWN'

        for dir in self.directions:
            parts = length_token.split()
            mesh_settings = mesh_settings + '\n' +  parts[0] + dir + ' ' + parts[1]

        mesh_settings = mesh_settings + '\n*DEPTH 1 1 1 12000.'

        self.write_content(mesh_settings)

    def write_effective_porosity(self):

        porosity_token = '*POR *ALL'
        
        # self.porosity_data = self.original_dataset[index_start:index_end]


    def write_effective_permeability(self):
        pass

    def write_content(self, content):

        with open('coarse_model.dat', 'a') as new:
            new.writelines(content)
