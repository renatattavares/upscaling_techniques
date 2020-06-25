import yaml
import copy
import numpy as np
from scipy.sparse import lil_matrix
from upscaling_procedures.local.local_upscaling import LocalUpscaling

class ChessBoard(LocalUpscaling):

    def set_simulation_variables(self):

        self.mesh.permeability[:] = np.array([7,7,7])
        self.mesh.permeability[np.array([0,1,4,5,16,17,20,21,10,11,26,27,31,30,14,15,45,44,60,61,56,57,40,41,50,51,34,35,54,55,39,38])] = np.array([3,3,3])
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
        correct_volumes_group_1 = np.array([0,1,2,3,16,17,18,19,32,33,34,35,48,49,50,51])
        correct_volumes_group_2 = np.array([12,13,14,15,28,29,30,31,44,45,46,47,60,61,62,63])
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
