import numpy as np
from scipy.sparse import lil_matrix

class ExtendedAssembly():

    def assembly_extended_local_problem(self, coarse_volume):

        volumes_global_ids, volumes_local_ids = self.volumes_extended_local_problem(coarse_volume)
        faces_global_ids, boundary_faces_global_ids, internal_faces_global_ids = self.faces_extended_local_problem(coarse_volume)
        self.equivalent_permeability = self.equivalent_permeability_extended_local_problem(faces_global_ids, internal_faces_global_ids)
        adjacent_volumes = self.mesh.faces.bridge_adjacencies(internal_faces_global_ids, 2, 3)
        adjacent_volumes_flatten = adjacent_volumes.flatten()
        global_ids_list = list(volumes_global_ids)

        # Rewrite adjacent_volumes array using local ids
        i = 0
        adjacent_volumes_local_ids = np.zeros(np.shape(adjacent_volumes), dtype = int)

        for duo in adjacent_volumes:
            new_duo = []
            for volume in duo:
                index = global_ids_list.index(volume)
                new_duo.append(index)
            adjacent_volumes_local_ids[i] = np.array(new_duo)
            i += 1

        neighbors_centers = np.reshape(self.mesh.volumes.center[adjacent_volumes_flatten], newshape = (len(internal_faces_global_ids), 6))
        neighbors_centers[:, 3:6] = neighbors_centers[:, 3:6]*(-1)
        centers_distance = np.linalg.norm((neighbors_centers[:, 0:3] + neighbors_centers[:, 3:6]), axis = 1); # Calculates the distante between the face's neighbors centers
        transmissibility = lil_matrix((int(len(volumes_global_ids))), (int(len(volumes_global_ids))), dtype = float)

        global_ids_list = list(faces_global_ids)
        internal_faces_local_ids = np.zeros(np.shape(internal_faces_global_ids), dtype = int)
        i = 0

        for id in internal_faces_global_ids:
            index = global_ids_list.index(id)
            internal_faces_local_ids[i] = index
            i += 1

        self.internal_faces_local_ids = internal_faces_local_ids
        self.internal_faces_global_ids = internal_faces_global_ids
        self.faces_global_ids = faces_global_ids
        self.adjacent_volumes_flatten = adjacent_volumes_flatten
        self.adjacent_volumes_local_ids = adjacent_volumes_local_ids

        # for j in range(len(internal_faces_local_ids)):
        #     id1 = int(adjacent_volumes_local_ids[j,0]) # ID of the first neighbor from the face
        #     id2 = int(adjacent_volumes_local_ids[j,1]) # ID of the second neighbor from the face
        #     transmissibility[id1,id2] += equivalent_permeability[internal_faces_local_ids[j]]/centers_distance[j]
        #     transmissibility[id2,id1] += equivalent_permeability[internal_faces_local_ids[j]]/centers_distance[j]

        lil_matrix.setdiag(transmissibility,(-1)*transmissibility.sum(axis = 1))

        return transmissibility

    def equivalent_permeability_extended_local_problem(self, faces_global_ids, internal_faces_global_ids):
        """
        Calculates the equivalent permeability of each internal face
        """
        equivalent_permeability = np.arange(len(faces_global_ids), dtype = float)

        for direction in self.direction_string:
            direction_number = self.directions_numbers.get(direction)
            direction_array = self.directions_dictionary.get(direction)
            parallel_direction = self.mesh.parallel_direction[internal_faces_global_ids]
            index_faces_direction = np.isin(parallel_direction, direction_number)
            index_faces_direction = np.where(index_faces_direction == True)[0]
            correct_faces = internal_faces_global_ids[index_faces_direction]
            adjacent_volumes = self.mesh.faces.bridge_adjacencies(correct_faces, 2, 3)
            permeability = self.mesh.permeability[adjacent_volumes.flatten()]
            permeability_direction = permeability[:,direction_number]
            permeability_direction = np.reshape(permeability_direction, newshape = (len(adjacent_volumes), 2))
            multiplication = np.multiply(permeability_direction[:,0], permeability_direction[:,1])
            sum = permeability_direction[:,0] + permeability_direction[:,1]
            global_ids_list = list(faces_global_ids)
            local_correct_faces = np.zeros(np.shape(correct_faces), dtype = int)
            i = 0

            for id in correct_faces:
                index = global_ids_list.index(id)
                local_correct_faces[i] = index
                i += 1

            equivalent_permeability[local_correct_faces] = 2*multiplication/sum

        return equivalent_permeability
