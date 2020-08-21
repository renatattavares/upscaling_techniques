"""
Module to rewrite specific keywords from IMEX data set after collecting upscaling informations of permeability, porosity and primal mesh
"""
import numpy as np
from imex_integration.refinement import MeshRefinement

class DatasetWriter(MeshRefinement):

    def __init__(self, original_dataset, coarsening, number_elements, length_elements, effective_porosity, effective_permeability, info_refined_volumes):

        print('\n##### Generating dataset file #####')

        np.set_printoptions(floatmode = 'fixed')
        np.set_printoptions(suppress = 'True')
        self.line = '\n'
        self.separator = '********************************************************************************'
        coarse_length = (number_elements/coarsening)*length_elements

        with open(original_dataset, 'r') as old:
            original_dataset = old.readlines()

        #self.get_default_settings(original_dataset, info_refined_volumes)

        with open('imex_datasets/coarse_model.dat', 'w') as new:
            pass

        print('Writing dataset heading...')
        #self.write_heading(original_dataset)
        print('Writing dataset mesh settings...')
        #self.write_mesh_settings(coarsening, coarse_length)
        #self.write_content(self.cpor)
        #self.write_content(self.prpor)

        if info_refined_volumes is not None:
            print('Writing refinement settings...')
            self.write_refine_card(length_elements, coarse_length, coarsening)

        print('Writing dataset effective porosity...')
        #self.write_effective_porosity(effective_porosity)
        print('Writing rock compressibility and reference pressure...')
        print('Writing dataset effective permeability...')
        #self.write_effective_permeability(effective_permeability)

        if info_refined_volumes is not None:
            self.write_refined_volumes_perm(info_refined_volumes, coarsening)
            self.write_refined_volumes_por(info_refined_volumes, coarsening)

        print('Writting model settings...')
        #self.write_content(self.model)

    def get_default_settings(self, original_dataset, info_refined_volumes):

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
                if info_refined_volumes is not None:
                    self.model = original_dataset[index:index+73]
                else:
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

    def write_refine_card(self, length_elements, coarse_length, coarsening):

        refine_token = '*REFINE'
        range_token = '*RANGE'
        self.upscaling_ratio = (coarse_length/length_elements).astype(int)
        refine_card = []

        refine_card.append(refine_token + ' ' + str(self.upscaling_ratio[0]) + ' ' + str(self.upscaling_ratio[1]) + ' ' + str(self.upscaling_ratio[2]))
        refine_card.append(self.line)
        refine_card.append(range_token + ' ' + '1 1 1:' + str(coarsening[2]))
        refine_card.append(self.line)
        refine_card.append(range_token + ' ' + '1 ' + str(coarsening[1]) + ' 1:' + str(coarsening[2]))
        refine_card.append(self.line)
        refine_card.append(range_token + ' ' + str(coarsening[0]) + ' 1 1:' + str(coarsening[2]))
        refine_card.append(self.line)
        refine_card.append(range_token + ' ' + str(coarsening[0]) + ' ' + str(coarsening[1]) + ' ' + '1:' + str(coarsening[2]))

        self.write_content(refine_card)

    def write_refined_volumes_perm(self, info_refined_volumes, coarsening):

        perm_direction_token = np.array(['*PERMI *RG', '*PERMJ *RG', '*PERMK *RG'])
        perm_refined_volumes = info_refined_volumes[0].flatten()
        perm_refined_volumes = np.reshape(perm_refined_volumes, newshape = (int(len(perm_refined_volumes)/3),3))
        direction = np.array([0,1,2])
        well1 = np.array([1,1,1], dtype = int)
        well2 = np.array([1,coarsening[1],1], dtype = int)
        well3 = np.array([coarsening[0],1,1], dtype = int)
        well4 = np.array([coarsening[0],coarsening[1],1], dtype = int)

        wells1 = []
        wells2 = []
        wells3 = []
        wells4 = []

        for i in range(coarsening[2]):
            wells1.append(well1 + np.array([0,0,i]))
            wells2.append(well2 + np.array([0,0,i]))
            wells3.append(well3 + np.array([0,0,i]))
            wells4.append(well4 + np.array([0,0,i]))

        perfs = np.array([wells1, wells2, wells3, wells4]).flatten()
        perfs = np.reshape(perfs, newshape = (int(len(perfs)/3),3))
        perm_card = []
        volumes_number = self.upscaling_ratio[0]*self.upscaling_ratio[1]*self.upscaling_ratio[2]

        for string, dir in zip(perm_direction_token, direction):
            perm = perm_refined_volumes[:,dir]
            i = 0
            for perf in perfs:
                p = perm[(0+(i*volumes_number)):(volumes_number+(i*volumes_number))]
                perm_card.append(string + ' ' + np.array2string(perf).replace('[', '').replace(']', '') + ' ' + '*ALL ' + np.array2string(p).replace('[', '').replace(']', ''))
                perm_card.append(self.line)
                i += 1
            perm_card.append(self.line)

        self.write_content(perm_card)

    def write_refined_volumes_por(self, info_refined_volumes, coarsening):
        por_token = '*POR *RG'
        por_refined_volumes = info_refined_volumes[1].flatten()
        well1 = np.array([1,1,1], dtype = int)
        well2 = np.array([1,coarsening[1],1], dtype = int)
        well3 = np.array([coarsening[0],1,1], dtype = int)
        well4 = np.array([coarsening[0],coarsening[1],1], dtype = int)

        wells1 = []
        wells2 = []
        wells3 = []
        wells4 = []

        for i in range(coarsening[2]):
            wells1.append(well1 + np.array([0,0,i]))
            wells2.append(well2 + np.array([0,0,i]))
            wells3.append(well3 + np.array([0,0,i]))
            wells4.append(well4 + np.array([0,0,i]))

        perfs = np.array([wells1, wells2, wells3, wells4]).flatten()
        perfs = np.reshape(perfs, newshape = (int(len(perfs)/3),3))
        por_card = []
        volumes_number = self.upscaling_ratio[0]*self.upscaling_ratio[1]*self.upscaling_ratio[2]

        i = 0
        for perf in perfs:
            p = por_refined_volumes[(0+(i*volumes_number)):(volumes_number+(i*volumes_number))]
            por_card.append(por_token + ' ' + np.array2string(perf).replace('[', '').replace(']', '') + ' ' + '*ALL ' + np.array2string(p).replace('[', '').replace(']', ''))
            por_card.append(self.line)
            i += 1

        self.write_content(por_card)

    def write_content(self, content):

        with open('imex_datasets/coarse_model.dat', 'a') as new:
            new.writelines(content)
            new.writelines(self.line)
