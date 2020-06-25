import yaml
import copy
import numpy as np
from scipy.sparse import lil_matrix
from upscaling_procedures.local.local_upscaling import LocalUpscaling

class Parallel(LocalUpscaling):

    def set_simulation_variables(self):

        self.mesh.permeability[np.array([2,3,6,7,10,11,14,15])] = np.array([3,3,3])
        self.mesh.permeability[np.array([0,1,4,5,8,9,12,13])] = np.array([7,7,7])
        self.mesh.porosity[:] = 1

        mesh_info_file = 'mesh_info.yml'

        with open(mesh_info_file, 'r') as file:
            data = yaml.safe_load(file)

        number_elements = data['Elements']
        lenght_elements = data['Length']
        self.number_elements = np.array([], dtype = int)
        self.length_elements = np.array([], dtype = int)

        for number, lenght in zip(number_elements.values(), lenght_elements.values()):
            self.number_elements = np.append(self.number_elements, number)
            self.length_elements = np.append(self.length_elements, lenght)

    def set_coordinate_system(self):
        self.x = np.array([1,0,0])
        self.y = np.array([0,1,0])
        self.z = np.array([0,0,1])
        self.direction_string = np.array(['z'])
        self.directions_dictionary = {
            'x': self.x,
            'y': self.y,
            'z': self.z
            }
        self.directions_numbers = {
            'x': 0,
            'y': 1,
            'z': 2
            }

    def fixed_constant_pressure(self, general_transmissibility):
        """
        Function to apply fixed constant pressure boundary condition, it returns a transmissibility and a source/sink matrix modified.
        """
        #print('Setting boundary conditions of local problem {}'.format(self.coarse_volume))
        correct_volumes_group_1 = np.array([0,1,2,3])
        correct_volumes_group_2 = np.array([12,13,14,15])
        transmissibility = copy.deepcopy(general_transmissibility)
        volumes_group_1 = correct_volumes_group_1
        volumes_group_2 = correct_volumes_group_2
        transmissibility[volumes_group_1] = 0
        transmissibility[volumes_group_2] = 0
        transmissibility[volumes_group_1, volumes_group_1] = 1
        transmissibility[volumes_group_2, volumes_group_2] = 1
        source = lil_matrix((int(self.number_volumes_local_problem), 1), dtype = 'float')
        source[volumes_group_1] = 1
        source[volumes_group_2] = 0

        #print('Fixed constant pressure boundary condition applied')

        return transmissibility, source, correct_volumes_group_1
