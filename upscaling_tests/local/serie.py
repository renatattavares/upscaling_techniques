import time
import yaml
import copy
import numpy as np
import multiprocessing as mp
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from upscaling_procedures.local.local_problems import LocalProblems
from impress.preprocessor.meshHandle.configTools.configClass import coarseningInit as coarse_config

class Serie(LocalProblems):

    def __init__(self, mesh_file = None, dataset = None):

        print('\n##### Treatment of local problems #####')

        # Preprocessing mesh with IMPRESS
        self.mesh_file = mesh_file
        self.preprocess_mesh()

        self.mesh.permeability[:] = np.array([1,1,1])
        self.mesh.porosity[:] = 1
        self.number_coarse_volumes = 1
        self.coarse_volume = 0
        self.number_volumes_local_problem = 12
        self.direction_string = np.array(['z'])
        self.x = np.array([1,0,0])
        self.y = np.array([0,1,0])
        self.z = np.array([0,0,1])
        self.directions_dictionary = {
            'x': self.x,
            'y': self.y,
            'z': self.z
            }

        self.solve_local_problems()

    def assembly_local_problem(self):
        """
        Assembly of local problems to generate the transmissibility matrix
        """
        i = self.coarse_volume

        faces_local_ids = self.coarse.elements[i].faces.internal # Local IDs from the internal faces from a coarse volume
        equivalent_permeability = self.equivalent_permeability(i, faces_local_ids)
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

        return transmissibility

    def equivalent_permeability(self, coarse_volume, faces_local_ids):
        """
        Calculates the equivalent permeability of each internal face
        """
        column = np.arange(3)
        equivalent_permeability = np.array((self.coarse.elements[coarse_volume].faces.all), dtype = float)

        for i, j in zip(self.direction_string, column):

            direction = self.directions_dictionary.get(i)
            global_faces_ids = self.coarse.elements[coarse_volume].faces.global_id[faces_local_ids]
            parallel_direction = self.mesh.parallel_direction[global_faces_ids]
            index_faces_direction = np.isin(parallel_direction, j)
            index_faces_direction = np.where(index_faces_direction == True)[0]
            correct_faces = faces_local_ids[index_faces_direction]
            adjacent_volumes = self.coarse.elements[coarse_volume].faces.bridge_adjacencies(correct_faces, 2, 3)
            global_id_adjacent_volumes = self.coarse.elements[coarse_volume].volumes.global_id[adjacent_volumes.flatten()]
            permeability = self.mesh.permeability[global_id_adjacent_volumes]
            permeability_direction = permeability[:,j]
            permeability_direction = np.reshape(permeability_direction, newshape = (len(adjacent_volumes), 2))
            multiplication = np.multiply(permeability_direction[:,0], permeability_direction[:,1])
            sum = permeability_direction[:,0] + permeability_direction[:,1]
            equivalent_permeability[correct_faces] = 2*multiplication/sum

        return equivalent_permeability

    def solve_local_problems(self):
        self.pressure_x = []
        self.pressure_y = []
        self.pressure_z = []

        for i in range(self.number_coarse_volumes):
            print('\nSolving local problem {} '.format(i))
            self.coarse_volume = i

            for j in self.direction_string:
                general_transmissibility = self.assembly_local_problem()
                print('{} direction'.format(j))
                transmissibility, source = self.fixed_constant_pressure(general_transmissibility)
                pressure = self.solver(transmissibility, source)

                if j == 'x':
                    self.pressure_x.append(pressure)
                elif j == 'y':
                    self.pressure_y.append(pressure)
                elif j == 'z':
                    self.pressure_z.append(pressure)

                print('Pressure stored')

    def fixed_constant_pressure(self, general_transmissibility):
        volumes_group_1 = np.array([0,1])
        volumes_group_2 = np.array([10,11])

        transmissibility = general_transmissibility
        transmissibility[volumes_group_1] = 0
        transmissibility[volumes_group_2] = 0
        transmissibility[volumes_group_1, volumes_group_1] = 1
        transmissibility[volumes_group_2, volumes_group_2] = 1
        source = lil_matrix((int(self.number_volumes_local_problem), 1), dtype = 'float')
        source[volumes_group_1] = 1
        source[volumes_group_2] = 0

        return transmissibility, source
