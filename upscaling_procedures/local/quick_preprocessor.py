import yaml
import numpy as np

class QuickPreprocessor:

    def get_mesh_informations(self, coarse_config, mesh_info_file = 'mesh_info.yml'):
        """
        Access coarsening informations given to IMPRESS.
        """
        tree = coarse_config.tree['Simple']
        self.nx = tree['nx']
        self.ny = tree['ny']
        self.nz = tree['nz']

        with open(mesh_info_file, 'r') as file:
            data = yaml.safe_load(file)

        self.number_elements_x_direction = data['x']
        self.number_elements_y_direction = data['y']
        self.number_elements_z_direction = data['z']
        self.number_coarse_volumes = len(self.coarse.elements) # Number of volumes from the coarse mesh
        self.number_volumes_local_problem = len(self.mesh.volumes)/(self.nx*self.ny*self.nz) # Number of fine scale volumes inside a coarse volume

        print('\nMesh informations accessed')

    def identify_top_bottom_volumes(self):

        if np.array_equal(self.direction, self.x) == True:
            direction = 0
        elif np.array_equal(self.direction, self.y) == True:
            direction = 1
        elif np.array_equal(self.direction, self.z) == True:
            direction = 2

        i = self.coarse_volume
        correct_volumes_group_1 = np.zeros((1, self.number_faces_coarse_face), dtype = int)
        correct_volumes_group_2 = np.zeros((1, self.number_faces_coarse_face), dtype = int)

        boundary_faces = self.coarse.elements[i].faces.boundary # Local IDs of boundary faces of a coarse volume
        global_ids_faces = self.coarse.elements[i].faces.global_id[boundary_faces]
        parallel_direction = self.mesh.parallel_direction[global_ids_faces]
        index_faces_direction = np.isin(parallel_direction, direction)
        index_faces_direction = np.where(index_faces_direction == True)[0]
        correct_faces = boundary_faces[index_faces_direction]

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

        global_ids_faces = self.coarse.elements[i].faces.global_id[:]
        index_group_1 = np.isin(global_ids_faces, group_1)
        index_group_2 = np.isin(global_ids_faces, group_2)
        local_ids_group_1 = np.where(index_group_1 == True)[0]
        local_ids_group_2 = np.where(index_group_2 == True)[0]

        correct_volumes_group_1 = self.coarse.elements[i].faces.bridge_adjacencies(local_ids_group_1, 2, 3).flatten()
        correct_volumes_group_2 = self.coarse.elements[i].faces.bridge_adjacencies(local_ids_group_2, 2, 3).flatten()

        if np.array_equal(self.direction, self.x) is True:
            self.correct_volumes_x = np.append(self.correct_volumes_x, correct_volumes_group_1)
        elif np.array_equal(self.direction, self.y) is True:
            self.correct_volumes_y = np.append(self.correct_volumes_y, correct_volumes_group_1)
        elif np.array_equal(self.direction, self.z) is True:
            self.correct_volumes_z = np.append(self.correct_volumes_z, correct_volumes_group_1)

        return correct_volumes_group_1, correct_volumes_group_2

    def identify_side_volumes(self):
        pass

    def check_parallel_direction(self):
        """
        Identifies the direction that is parallel to the normal vector of each face
        """
        self.mesh.parallel_direction[:] = 4
        j = 0

        for i in self.direction_string:
            self.direction = self.directions_dictionary.get(i)
            global_ids_faces = np.arange(len(self.mesh.faces))
            normal_faces = self.mesh.faces.normal[global_ids_faces]
            direction_vector = np.ndarray(shape = np.shape(normal_faces), dtype = float)
            direction_vector = np.full_like(direction_vector, self.direction)

            # Cross product and norm
            cross_product = np.cross(normal_faces, direction_vector)
            norm_cross_product = np.linalg.norm(cross_product, axis = 1)

            # Verify which norms are one (if norm == 1, the face is parallel to the direction)
            correct_faces = np.isin(norm_cross_product, 0)
            index_correct_faces = np.where(correct_faces == True)[0]
            correct_faces = global_ids_faces[index_correct_faces]
            self.mesh.parallel_direction[correct_faces] = int(j)
            j+= 1
