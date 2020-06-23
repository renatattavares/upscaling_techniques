import numpy as np
from scipy.sparse import lil_matrix

class Assembly:

    def equivalent_permeability(self, faces_local_ids):
        """
        Calculates the equivalent permeability of each internal face
        """
        equivalent_permeability = np.array((self.coarse.elements[self.coarse_volume].faces.all), dtype = float)

        for direction in self.direction_string:
            direction_number = self.directions_numbers.get(direction)
            direction_array = self.directions_dictionary.get(direction)
            global_faces_ids = self.coarse.elements[self.coarse_volume].faces.global_id[faces_local_ids]
            parallel_direction = self.mesh.parallel_direction[global_faces_ids]
            index_faces_direction = np.isin(parallel_direction, direction_number)
            index_faces_direction = np.where(index_faces_direction == True)[0]
            correct_faces = faces_local_ids[index_faces_direction]
            adjacent_volumes = self.coarse.elements[self.coarse_volume].faces.bridge_adjacencies(correct_faces, 2, 3)
            global_id_adjacent_volumes = self.coarse.elements[self.coarse_volume].volumes.global_id[adjacent_volumes.flatten()]
            permeability = self.mesh.permeability[global_id_adjacent_volumes]
            permeability_direction = permeability[:,direction_number]
            permeability_direction = np.reshape(permeability_direction, newshape = (len(adjacent_volumes), 2))
            multiplication = np.multiply(permeability_direction[:,0], permeability_direction[:,1])
            sum = permeability_direction[:,0] + permeability_direction[:,1]
            equivalent_permeability[correct_faces] = 2*multiplication/sum

        return equivalent_permeability

    def assembly_local_problem(self):
        """
        Assembly of local problems to generate the transmissibility matrix
        """

        #print("Assembly of local problem {}".format(i))
        faces_local_ids = self.coarse.elements[self.coarse_volume].faces.internal # Local IDs from the internal faces from a coarse volume
        equivalent_permeability = self.equivalent_permeability(faces_local_ids)
        adjacent_volumes = self.coarse.elements[self.coarse_volume].faces.bridge_adjacencies(faces_local_ids, 2, 3) # Local IDs from both the neighbors from each of the internal faces
        adjacent_volumes_flatten = adjacent_volumes.flatten()
        neighbors_centers = np.reshape(self.coarse.elements[self.coarse_volume].volumes.center[adjacent_volumes_flatten], newshape = (len(faces_local_ids), 6))
        neighbors_centers[:, 3:6] = neighbors_centers[:, 3:6]*(-1)
        centers_distance = np.linalg.norm((neighbors_centers[:, 0:3] + neighbors_centers[:, 3:6]), axis = 1); # Calculates the distante between the face's neighbors centers
        transmissibility = lil_matrix((int(self.number_volumes_local_problem), int(self.number_volumes_local_problem)), dtype = float)

        for j in range(len(faces_local_ids)):
            id1 = int(adjacent_volumes[j,0]) # ID of the first neighbor from the face
            id2 = int(adjacent_volumes[j,1]) # ID of the second neighbor from the face
            transmissibility[id1,id2] += equivalent_permeability[faces_local_ids[j]]/centers_distance[j]
            transmissibility[id2,id1] += equivalent_permeability[faces_local_ids[j]]/centers_distance[j]

        lil_matrix.setdiag(transmissibility,(-1)*transmissibility.sum(axis = 1))

        return transmissibility
