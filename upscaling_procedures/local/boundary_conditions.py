import copy
import numpy as np
from scipy.sparse import lil_matrix
from upscaling_procedures.local.quick_preprocessor import QuickPreprocessor

class BoundaryConditions(QuickPreprocessor):

    def __init__(self, boundary_condition_type):
        self.boundary_condition_type = boundary_condition_type

    def set_boundary_conditions(self, general_transmissibility):
        """
        Indicates which function must be executed to set boundary conditions acording to the option informed by the user.
        """

        self.boundary_conditions_dictionary = {
            1: self.fixed_constant_pressure, # Fixed constant pressure
            2: self.fixed_linear_pressure,   # Fixed linear pressure
            3: self.periodic_pressure        # Periodic pressure
            }

        self.check_directions_and_coarse_face()

        general_transmissibility = general_transmissibility
        transmissibility, source, correct_volumes_group_1 = self.boundary_conditions_dictionary.get(self.boundary_condition_type, "\nprint('Invalid boundary condition')")(general_transmissibility) # Execute the correct boundary condition function

        return transmissibility, source, correct_volumes_group_1

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

        return transmissibility, source, correct_volumes_group_1

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
