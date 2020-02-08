import numpy as np
from scipy.sparse import lil_matrix

class Assembly:

    def equivalent_permeability(self, coarse_volume, faces_local_ids):
        """
        Calculates the equivalent permeability of each internal face
        """
        faces_normal = self.coarse.elements[coarse_volume].faces.normal[faces_local_ids]
        column = np.arange(3)
        equivalent_permeability = lil_matrix((len(self.coarse.elements[coarse_volume].faces.all), 1), dtype = float)

        for i, j in zip(self.direction_string, column):
            direction = self.directions_dictionary.get(i)
            cross_product = np.cross(direction, faces_normal)
            norm_cross_product = np.linalg.norm(cross_product, axis = 1)
            correct_faces = np.isin(norm_cross_product, 0)
            index_faces_direction = np.where(correct_faces == 1)[0]
            faces_direction = faces_local_ids[index_faces_direction]
            adjacent_volumes = self.coarse.elements[coarse_volume].faces.bridge_adjacencies(faces_direction, 2, 3)
            global_id_adjacent_volumes = self.coarse.elements[coarse_volume].volumes.global_id[adjacent_volumes]
            permeability_direction = self.mesh.permeability[global_id_adjacent_volumes]
            permeability_direction = permeability_direction[:,j]
            permeability_direction = np.reshape(permeability_direction, newshape = (len(adjacent_volumes), 2))
            p1 = permeability_direction[:,0]
            p2 = permeability_direction[:,1]
            permeability_multiplication = np.multiply(p1, p2)
            permeability_sum = p1 + p2
            result = 2*permeability_multiplication/permeability_sum
            equivalent_permeability[faces_direction] = np.reshape(result, newshape = (len(faces_direction), 1))

        return equivalent_permeability


    def assembly_local_problem(self):
        """
        Assembly of local problems to generate the transmissibility matrix
        """
        self.transmissibilities = np.array([])

        for i in range(self.number_coarse_volumes):
            print("Assembly of local problem {}".format(i))
            faces_local_ids = self.coarse.elements[i].faces.internal # Local IDs from the internal faces from a coarse volume
            equivalent_permeability = self.equivalent_permeability(i, faces_local_ids)
            equivalent_permeability = lil_matrix.toarray(equivalent_permeability)
            adjacent_volumes = self.coarse.elements[i].faces.bridge_adjacencies(faces_local_ids, 2, 3) # Local IDs from both the neighbors from each of the internal faces
            adjacent_volumes_flatten = adjacent_volumes.flatten()
            neighbors_centers = np.reshape(self.coarse.elements[i].volumes.center[adjacent_volumes_flatten], newshape = (len(faces_local_ids), 6))
            neighbors_centers[:, 3:6] = neighbors_centers[:, 3:6]*(-1)
            centers_distance = np.linalg.norm((neighbors_centers[:, 0:3] + neighbors_centers[:, 3:6]), axis = 1); # Calculates the distante between the face's neighbors centers
            transmissibility = lil_matrix((int(self.number_volumes_local_problem), int(self.number_volumes_local_problem)), dtype = float)

            for j in range(len(faces_local_ids)):
                id1 = int(adjacent_volumes[j,0]) # ID of the first neighbor from the face
                id2 = int(adjacent_volumes[j,1]) # ID of the second neighbor from the face
                transmissibility[id1,id2] += equivalent_permeability[faces_local_ids[j]]/centers_distance[j]
                transmissibility[id2,id1] += equivalent_permeability[faces_local_ids[j]]/centers_distance[j]

            lil_matrix.setdiag(transmissibility,(-1)*transmissibility.sum(axis = 1))
            self.transmissibilities = np.append(self.transmissibilities, transmissibility) # Store transmissibility
