import numpy as np
from scipy.sparse import lil_matrix

class BoundaryConditionsLocalProblems:

    def __init__(self, data, boundary_condition_type):
        self.data = data
        self.boundary_condition_type = boundary_condition_type

    def set_boundary_conditions(self, direction):
        """
        Indicates which function must be executed to set boundary conditions acording to the option informed by the user.
        """
        
        self.data.boundary_conditions_dictionary = {
            1: self.fixed_constant_pressure, # Fixed constant pressure
            2: self.fixed_linear_pressure,   # Fixed linear pressure
            3: self.periodic_pressure        # Periodic pressure
            }

        self.data.direction = self.data.directions_dictionary.get(direction) # Get the direction given
        if np.array_equal(self.data.direction, self.data.x) == True:
            self.perpendicular_direction_1 = self.data.y
            self.perpendicular_direction_2 = self.data.z
            self.number_faces_coarse_face = int((self.data.number_elements_y_direction/self.data.ny)*(self.data.number_elements_z_direction/self.data.nz))

        elif np.array_equal(self.data.direction, self.data.y) == True:
            self.perpendicular_direction_1 = self.data.x
            self.perpendicular_direction_2 = self.data.z
            self.number_faces_coarse_face = int((self.data.number_elements_x_direction/self.data.nx)*(self.data.number_elements_z_direction/self.data.nz))

        elif np.array_equal(self.data.direction, self.data.z) == True:
            self.perpendicular_direction_1 = self.data.x
            self.perpendicular_direction_2 = self.data.y
            self.number_faces_coarse_face = int((self.data.number_elements_x_direction/self.data.ny)*(self.data.number_elements_y_direction/self.data.ny))

        self.data.boundary_conditions_dictionary.get(self.boundary_condition_type, "\nprint('Invalid boundary condition')")() # Execute the correct boundary condition function

    def fixed_constant_pressure(self):
        """
        Function to apply fixed constant pressure boundary condition, it returns a transmissibility and a source/sink matrix modified.
        """
        correct_volumes_group_1, correct_volumes_group_2 = self.identify_top_bottom_volumes()

        for i in range(self.data.number_coarse_volumes):
            volumes_group_1 = correct_volumes_group_1[i]
            volumes_group_2 = correct_volumes_group_2[i]
            self.data.coarse.elements[i].transmissibility[volumes_group_1] = 0
            self.data.coarse.elements[i].transmissibility[volumes_group_2] = 0
            self.data.coarse.elements[i].transmissibility[volumes_group_1, volumes_group_1] = 1
            self.data.coarse.elements[i].transmissibility[volumes_group_2, volumes_group_2] = 1
            self.data.coarse.elements[i].source = lil_matrix((int(self.data.number_volumes_local_problem), 1), dtype = 'float')
            self.data.coarse.elements[i].source[volumes_group_1] = self.data.pressure_gradient
            self.data.coarse.elements[i].source[volumes_group_2] = 0

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

    def identify_top_bottom_volumes(self):

        correct_volumes_group_1 = np.zeros((self.data.number_coarse_volumes, self.number_faces_coarse_face), dtype = int)
        correct_volumes_group_2 = np.zeros((self.data.number_coarse_volumes, self.number_faces_coarse_face), dtype = int)


        for i in range(self.data.number_coarse_volumes):
            boundary_faces = self.data.coarse.elements[i].faces.boundary # Local IDs of boundary faces of a coarse volume
            normal_boundary_faces = self.data.coarse.elements[i].faces.normal[boundary_faces] # Normal vector of boundary faces of a coarse volumes
            direction_vector = np.ndarray(shape = np.shape(normal_boundary_faces), dtype = float)
            direction_vector = np.full_like(direction_vector, self.data.direction) # Vectorization

            # Cross product and norm
            cross_product = np.cross(normal_boundary_faces, direction_vector)
            norm_cross_product = np.linalg.norm(cross_product, axis = 1)

            # Verify which norms are zero (if norm == 0, the face is perpendicular to the direction)
            correct_faces = np.isin(norm_cross_product, 0)
            index_correct_faces = np.where(correct_faces == True)[0]
            correct_faces = boundary_faces[index_correct_faces]

            # Separate faces in two groups
            global_ids_correct_faces = self.data.coarse.elements[i].faces.global_id[correct_faces]
            interface_coarse_face_id = self.data.coarse.iface_neighbors(i)[1]

            for j in range(len(interface_coarse_face_id)):
                interface_faces = self.data.coarse.interfaces_faces[int(interface_coarse_face_id[j])]
                verify = np.any(np.isin(interface_faces, global_ids_correct_faces[0]))
                if verify == True:
                    index = np.isin(interface_faces, global_ids_correct_faces)
                    group_1 = interface_faces[index]
                    break

            index = np.isin(global_ids_correct_faces, group_1, assume_unique = False, invert = True)
            group_2 = global_ids_correct_faces[index]

            global_ids_faces = self.data.coarse.elements[i].faces.father_id[:]
            index_group_1 = np.isin(global_ids_faces, group_1)
            index_group_2 = np.isin(global_ids_faces, group_2)
            local_ids_group_1 = np.where(index_group_1 == True)[0]
            local_ids_group_2 = np.where(index_group_2 == True)[0]

            adjacent_volumes_group_1 = self.data.coarse.elements[i].faces.bridge_adjacencies(local_ids_group_1, 2, 3)
            adjacent_volumes_group_1 = np.reshape(adjacent_volumes_group_1, newshape = (1,25))
            correct_volumes_group_1[i] = adjacent_volumes_group_1
            adjacent_volumes_group_2 = self.data.coarse.elements[i].faces.bridge_adjacencies(local_ids_group_2, 2, 3)
            adjacent_volumes_group_2 = np.reshape(adjacent_volumes_group_2, newshape = (1,25))
            correct_volumes_group_2[i] = adjacent_volumes_group_2

        return correct_volumes_group_1, correct_volumes_group_2

    def identify_side_volumes(self):
        pass
