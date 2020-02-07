import yaml
import numpy as np

class QuickPreprocessor:

    def get_mesh_informations(self, mesh_info_file = 'mesh_info.yml'):
        """
        Access coarsening informations given to IMPRESS.
        """
        self.tree = self.coarse_config.tree['Simple']
        self.nx = self.tree['nx']
        self.ny = self.tree['ny']
        self.nz = self.tree['nz']

        with open(mesh_info_file, 'r') as file:
            data = yaml.safe_load(file)

        self.number_elements_x_direction = data['x']
        self.number_elements_y_direction = data['y']
        self.number_elements_z_direction = data['z']

        print('\nMesh informations accessed')

    def identify_top_bottom_volumes(self):

        self.correct_volumes_group_1 = np.zeros((self.number_coarse_volumes, self.number_faces_coarse_face), dtype = int)
        self.correct_volumes_group_2 = np.zeros((self.number_coarse_volumes, self.number_faces_coarse_face), dtype = int)

        for i in range(self.number_coarse_volumes):
            boundary_faces = self.coarse.elements[i].faces.boundary # Local IDs of boundary faces of a coarse volume
            normal_boundary_faces = self.coarse.elements[i].faces.normal[boundary_faces] # Normal vector of boundary faces of a coarse volumes
            direction_vector = np.ndarray(shape = np.shape(normal_boundary_faces), dtype = float)
            direction_vector = np.full_like(direction_vector, self.direction) # Vectorization

            # Cross product and norm
            cross_product = np.cross(normal_boundary_faces, direction_vector)
            norm_cross_product = np.linalg.norm(cross_product, axis = 1)

            # Verify which norms are zero (if norm == 0, the face is perpendicular to the direction)
            correct_faces = np.isin(norm_cross_product, 0)
            index_correct_faces = np.where(correct_faces == True)[0]
            correct_faces = boundary_faces[index_correct_faces]

            # Separate faces in two groups
            global_ids_correct_faces = self.coarse.elements[i].faces.global_id[correct_faces]
            interface_coarse_face_id = self.coarse.iface_neighbors(i)[1]

            for j in range(len(interface_coarse_face_id)):
                interface_faces = self.coarse.interfaces_faces[int(interface_coarse_face_id[j])]
                verify = np.any(np.isin(interface_faces, global_ids_correct_faces[0]))
                if verify == True:
                    index = np.isin(interface_faces, global_ids_correct_faces)
                    group_1 = interface_faces[index]
                    break

            index = np.isin(global_ids_correct_faces, group_1, assume_unique = False, invert = True)
            group_2 = global_ids_correct_faces[index]

            global_ids_faces = self.coarse.elements[i].faces.father_id[:]
            index_group_1 = np.isin(global_ids_faces, group_1)
            index_group_2 = np.isin(global_ids_faces, group_2)
            local_ids_group_1 = np.where(index_group_1 == True)[0]
            local_ids_group_2 = np.where(index_group_2 == True)[0]

            adjacent_volumes_group_1 = self.coarse.elements[i].faces.bridge_adjacencies(local_ids_group_1, 2, 3).flatten()
            #adjacent_volumes_group_1 = np.reshape(adjacent_volumes_group_1, newshape = (1,self.number_faces_coarse_face))
            self.correct_volumes_group_1[i] = adjacent_volumes_group_1
            adjacent_volumes_group_2 = self.coarse.elements[i].faces.bridge_adjacencies(local_ids_group_2, 2, 3).flatten()
            #adjacent_volumes_group_2 = np.reshape(adjacent_volumes_group_2, newshape = (1,self.number_faces_coarse_face))
            self.correct_volumes_group_2[i] = adjacent_volumes_group_2

        if np.array_equal(self.direction, self.x) is True:
            self.correct_volumes_x = self.correct_volumes_group_1
        elif np.array_equal(self.direction, self.y) is True:
            self.correct_volumes_y = self.correct_volumes_group_1
        elif np.array_equal(self.direction, self.z) is True:
            self.correct_volumes_z = self.correct_volumes_group_1

    def identify_side_volumes(self):
        pass
