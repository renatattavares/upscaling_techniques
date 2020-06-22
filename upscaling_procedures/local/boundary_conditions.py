import copy
import numpy as np
from scipy.sparse import lil_matrix
from upscaling_procedures.local.quick_preprocessor import QuickPreprocessor

class BoundaryConditions(QuickPreprocessor):

    def __init__(self, boundary_condition_type):
        self.boundary_condition_type = boundary_condition_type

    def set_boundary_conditions(self, direction, general_transmissibility):
        """
        Indicates which function must be executed to set boundary conditions acording to the option informed by the user.
        """

        self.boundary_conditions_dictionary = {
            1: self.fixed_constant_pressure, # Fixed constant pressure
            2: self.fixed_linear_pressure,   # Fixed linear pressure
            3: self.periodic_pressure        # Periodic pressure
            }

        self.direction = self.directions_dictionary.get(direction) # Get the direction given

        if np.array_equal(self.direction, self.x) is True:
            self.perpendicular_direction_1 = self.y
            self.perpendicular_direction_2 = self.z
            self.number_faces_coarse_face = int((self.number_elements_y_direction/self.ny)*(self.number_elements_z_direction/self.nz))

        elif np.array_equal(self.direction, self.y) is True:
            self.perpendicular_direction_1 = self.x
            self.perpendicular_direction_2 = self.z
            self.number_faces_coarse_face = int((self.number_elements_x_direction/self.nx)*(self.number_elements_z_direction/self.nz))

        elif np.array_equal(self.direction, self.z) is True:
            self.perpendicular_direction_1 = self.x
            self.perpendicular_direction_2 = self.y
            self.number_faces_coarse_face = int((self.number_elements_x_direction/self.ny)*(self.number_elements_y_direction/self.ny))

        general_transmissibility = general_transmissibility
        transmissibility, source = self.boundary_conditions_dictionary.get(self.boundary_condition_type, "\nprint('Invalid boundary condition')")(general_transmissibility) # Execute the correct boundary condition function

        return transmissibility, source

    def fixed_constant_pressure(self, general_transmissibility):
        """
        Function to apply fixed constant pressure boundary condition, it returns a transmissibility and a source/sink matrix modified.
        """
        #print('Setting boundary conditions of local problem {}'.format(self.coarse_volume))
        correct_volumes_group_1, correct_volumes_group_2 = self.identify_top_bottom_volumes()
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

        return transmissibility, source

    # Depending on the identify_side_volumes function
    def fixed_linear_pressure(self):
        """
        Function to apply fixed linear pressure boundary condition, it returns a transmissibility and a source/sink matrix modified.
        """

        print('\nFixed linear pressure boundary condition applied')

    def periodic_pressure(self):
        """
        Function to apply periodic pressure boundary condition, it returns a transmissibility and a source/sink matrix modified.
        """

        print('\nPeriodic pressure boundary condition applied')
