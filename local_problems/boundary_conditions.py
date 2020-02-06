import numpy as np
from scipy.sparse import lil_matrix
from local_problems.quick_preprocessor import QuickPreprocessor

class BoundaryConditions(QuickPreprocessor):

    def __init__(self, boundary_condition_type):
        self.boundary_condition_type = boundary_condition_type

    def set_boundary_conditions(self, direction):
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

        self.boundary_conditions_dictionary.get(self.boundary_condition_type, "\nprint('Invalid boundary condition')")() # Execute the correct boundary condition function

    def fixed_constant_pressure(self):
        """
        Function to apply fixed constant pressure boundary condition, it returns a transmissibility and a source/sink matrix modified.
        """
        self.identify_top_bottom_volumes()

        for i in range(self.number_coarse_volumes):
            volumes_group_1 = self.correct_volumes_group_1[i]
            volumes_group_2 = self.correct_volumes_group_2[i]
            self.coarse.elements[i].transmissibility[volumes_group_1] = 0
            self.coarse.elements[i].transmissibility[volumes_group_2] = 0
            self.coarse.elements[i].transmissibility[volumes_group_1, volumes_group_1] = 1
            self.coarse.elements[i].transmissibility[volumes_group_2, volumes_group_2] = 1
            self.coarse.elements[i].source = lil_matrix((int(self.number_volumes_local_problem), 1), dtype = 'float')
            self.coarse.elements[i].source[volumes_group_1] = self.pressure_gradient
            self.coarse.elements[i].source[volumes_group_2] = 0

        print('\nFixed constant pressure boundary condition applied')

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
