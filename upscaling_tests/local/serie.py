import yaml
import numpy as np
from scipy.sparse import lil_matrix
from upscaling_procedures.local.local_upscaling import LocalUpscaling

class Serie(LocalUpscaling):

    def set_simulation_variables(self):

        self.mesh.permeability[np.array([0,1,2,3,4,5])] = np.array([7,7,7])
        self.mesh.permeability[np.array([6,7,8,9,10,11])] = np.array([3,3,3])
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

    def upscale_permeability_porosity(self):
        """
        It calculates the effective permeability of a coarse volume
        """

        self.effective_permeability = []
        self.effective_porosity = []
        volume = self.length_elements[0]*self.length_elements[1]*self.length_elements[2]

        for cv in range(self.number_coarse_volumes):
            print('Upscaling of local problem {}'.format(cv))
            self.coarse_volume = cv
            general_transmissibility = self.assembly_local_problem()
            effective_permeability = []
            total_volume = volume*self.number_volumes_local_problem
            global_ids_volumes = self.coarse.elements[cv].volumes.father_id[:]
            porosity = self.mesh.porosity[global_ids_volumes]
            effective_porosity = (volume*porosity).sum()/total_volume

            for direction in self.direction_string:
                #print('in {} direction'.format(direction))
                self.direction = direction
                transmissibility, source, local_wall = self.set_boundary_conditions(general_transmissibility)
                pressures = self.solver(transmissibility, source)

                direction_number = self.directions_numbers.get(direction)
                center_distance_walls = self.center_distance_walls[cv][direction_number]
                area = self.areas[direction_number]
                global_wall = self.coarse.elements[cv].volumes.global_id[local_wall]
                local_adj = self.identify_adjacent_volumes_to_wall(cv, local_wall)
                global_adj = self.coarse.elements[cv].volumes.global_id[local_adj]
                center_wall = self.coarse.elements[cv].volumes.center[local_wall]
                center_adj = self.coarse.elements[cv].volumes.center[local_adj]
                pressure_wall = pressures[local_wall]
                pressure_adj = pressures[local_adj]
                pw = self.mesh.permeability[global_wall]
                permeability_wall = pw[:, direction_number]
                pa = self.mesh.permeability[global_adj]
                permeability_adj = pa[:, direction_number]
                flow_rate = ((((2*np.multiply(permeability_wall,permeability_adj)/(permeability_wall+permeability_adj))*(pressure_wall-pressure_adj))*self.areas[direction_number])/np.linalg.norm(center_wall - center_adj, axis = 1)).sum()
                effective_permeability.append(center_distance_walls*flow_rate/(area*self.number_faces_coarse_face[direction_number]))

            self.effective_permeability.append(effective_permeability)
            self.effective_porosity.append(effective_porosity)

    def fixed_constant_pressure(self, general_transmissibility):
        """
        Function to apply fixed constant pressure boundary condition, it returns a transmissibility and a source/sink matrix modified.
        """
        #print('Setting boundary conditions of local problem {}'.format(self.coarse_volume))
        correct_volumes_group_1 = np.array([0,1])
        correct_volumes_group_2 = np.array([10,11])
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
