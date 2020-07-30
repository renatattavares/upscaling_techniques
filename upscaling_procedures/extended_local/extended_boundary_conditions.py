import copy
import numpy as np
from scipy.sparse import lil_matrix
from upscaling_procedures.local.boundary_conditions import BoundaryConditions

class ExtendedBoundaryConditions(BoundaryConditions):

    def set_boundary_conditions(self, general_transmissibility, coarse_volume):
        """
        Indicates which function must be executed to set boundary conditions acording to the option informed by the user.
        """

        self.boundary_conditions_dictionary = {
            1: self.fixed_constant_pressure, # Fixed constant pressure
            2: self.fixed_linear_pressure,   # Fixed linear pressure
            3: self.periodic_pressure        # Periodic pressure
            }

        coarse_volume = coarse_volume
        general_transmissibility = general_transmissibility
        transmissibility, source, correct_volumes_group_1 = self.boundary_conditions_dictionary.get(self.boundary_condition_type, "\nprint('Invalid boundary condition')")(general_transmissibility, coarse_volume) # Execute the correct boundary condition function

        return transmissibility, source, correct_volumes_group_1

    def fixed_constant_pressure(self, general_transmissibility, coarse_volume):
        """
        Function to apply fixed constant pressure boundary condition, it returns a transmissibility and a source/sink matrix modified.
        """
        #print('Setting boundary conditions of local problem {}'.format(self.coarse_volume))
        correct_volumes_group_1, correct_volumes_group_2 = self.identify_top_bottom_volumes(coarse_volume)
        volumes_global_ids, volumes_local_ids = self.volumes_extended_local_problem(coarse_volume)
        transmissibility = copy.deepcopy(general_transmissibility)
        volumes_group_1 = correct_volumes_group_1
        volumes_group_2 = correct_volumes_group_2
        transmissibility[volumes_group_1] = 0
        transmissibility[volumes_group_2] = 0
        transmissibility[volumes_group_1, volumes_group_1] = 1
        transmissibility[volumes_group_2, volumes_group_2] = 1
        source = lil_matrix((int(len(volumes_global_ids)), 1), dtype = 'float')
        source[volumes_group_1] = 1
        source[volumes_group_2] = 0

        #print('Fixed constant pressure boundary condition applied')

        return transmissibility, source, correct_volumes_group_1

    def identify_top_bottom_volumes(self, coarse_volume):

        direction = self.directions_dictionary.get(self.direction)
        direction_number = self.directions_numbers.get(self.direction)

        volumes_global_ids, volumes_local_ids = self.volumes_extended_local_problem(coarse_volume)
        volumes_centers = self.mesh.volumes.center[volumes_global_ids]
        coords = volumes_centers[:, direction_number]
        greatest = np.amax(coords)
        smallest = np.amin(coords)
        correct_volumes_group_1 = np.where(coords == smallest)[0]
        correct_volumes_group_2 = np.where(coords == greatest)[0]
        correct_volumes_group_1 = volumes_global_ids[correct_volumes_group_1]
        correct_volumes_group_2 = volumes_global_ids[correct_volumes_group_2]
        correct_volumes_group_1 = self.global_to_local_id(volumes_global_ids, correct_volumes_group_1)
        correct_volumes_group_2 = self.global_to_local_id(volumes_global_ids, correct_volumes_group_2)
        self.number_volumes_local_problem = len(volumes_global_ids)

        return correct_volumes_group_1, correct_volumes_group_2
