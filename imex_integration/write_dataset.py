"""
Module to rewrite specific keywords from IMEX data set after collecting upscaling informations of permeability, porosity and primal mesh
"""
import numpy as np

class DatasetWriter:

    def __init__(self, original_dataset, coarsening, length_elements, effective_porosity, effective_permeability):

        print('\n##### Generating dataset file #####')

        np.set_printoptions(floatmode = 'fixed')
        np.set_printoptions(suppress = 'True')
        self.line = '\n'
        self.separator = '********************************************************************************'

        with open(original_dataset, 'r') as old:
            original_dataset = old.readlines()

        self.get_default_settings(original_dataset)

        with open('imex_datasets/coarse_model.dat', 'w') as new:
            pass

        print('Writing dataset heading...')
        self.write_heading(original_dataset)
        print('Writing dataset mesh settings...')
        self.write_mesh_settings(coarsening, length_elements)
        print('Writing dataset effective porosity...')
        self.write_effective_porosity(effective_porosity)
        print('Writing rock compressibility and reference pressure...')
        self.write_content(self.cpor)
        self.write_content(self.prpor)
        print('Writing dataset effective permeability...')
        self.write_effective_permeability(effective_permeability)
        print('Writting model settings...')
        self.write_content(self.model)

    def get_default_settings(self, original_dataset):

        cpor = '*CPOR'
        prpor = '*PRPOR'
        model = '*MODEL *OILWATER'

        for line in original_dataset:
            if cpor in line:
                self.cpor = line
            if prpor in line:
                self.prpor = line
            if model in line:
                index = original_dataset.index(line)
                self.model = original_dataset[index:]

    def write_heading(self, original_dataset):

        heading_token = '**         Reservoir Description Section                                      **'

        for line in original_dataset:
            if heading_token in line:
                index = original_dataset.index(line)

        heading = original_dataset[:index]
        heading.append(heading_token)
        heading.append(self.line)
        heading.append(self.separator)
        self.write_content(heading)

    def write_mesh_settings(self, coarsening, length):

        grid_token = '*GRID *CART'
        length_token = ['*DI *CON', '*DJ *CON', '*DK *CON']
        mesh_settings = []
        values_line = ''

        for value in coarsening.astype(str):
            values_line = values_line + ' ' + value

        values_line = grid_token + values_line
        mesh_settings.append(values_line)
        mesh_settings.append(self.line)
        mesh_settings.append('*KDIR *DOWN')
        mesh_settings.append(self.line)

        for value, token in zip(length.astype(str), length_token):
            mesh_settings.append(token + ' ' + value)
            mesh_settings.append(self.line)

        mesh_settings.append('*DEPTH 1 1 1 12000.')
        self.write_content(mesh_settings)

    def write_effective_porosity(self, effective_porosity):

        porosity_token = '*POR *ALL'
        porosity_card = []
        porosity_card.append(porosity_token)
        porosity_card.append(self.line)
        ep = np.reshape(effective_porosity, newshape = (int(len(effective_porosity)/5), 5))

        for line in ep:
            str_line = np.array2string(line).replace('[', '').replace(']', '')
            porosity_card.append(str_line)
            porosity_card.append(self.line)

        self.write_content(porosity_card)

    def write_effective_permeability(self, effective_permeability):

        perm_direction_token = np.array(['*PERMI *ALL', '*PERMJ *ALL', '*PERMK *ALL'])
        direction = np.array([0,1,2])

        for string, dir in zip(perm_direction_token, direction):
            permeability_card = []
            permeability_card.append(string)
            permeability_card.append(self.line)
            ep = np.reshape(effective_permeability[:,dir], newshape = (int(len(effective_permeability[:,dir])/5), 5))

            for line in ep:
                str_line = np.array2string(line).replace('[', '').replace(']', '')
                permeability_card.append(str_line)
                permeability_card.append(self.line)

            self.write_content(permeability_card)

    def write_content(self, content):

        with open('imex_datasets/coarse_model.dat', 'a') as new:
            new.writelines(content)
            new.writelines(self.line)
