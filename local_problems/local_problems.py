"""
Module for treatment of local problems to apply local upscaling technique in structured tridimensional meshes
"""
import time
import numpy as np
#import xlsxwriter
import yaml
from scipy.sparse import lil_matrix
from local_problems.solver_local_problems import Solver
from local_problems.boundary_conditions import BoundaryConditions
from impress.preprocessor.meshHandle.configTools.configClass import coarseningInit as coarse_config
from impress.preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as impress

class LocalProblems(BoundaryConditions, Solver):

    def __init__(self, mesh_file = None, boundary_condition_type = None):

        print('\n##### Treatment of local problems #####')

        # Preprocessing mesh with IMPRESS
        print('\nPre-processing mesh with IMPRESS...')
        start = time.time()
        self.mesh = impress(mesh_file, dim = 3) # IMPRESS' object
        end = time.time()
        print("\nThis step lasted {0}s".format(end-start))

        # Setting variables
        print('\nAccessing coarsening informations from IMPRESS and setting important variables...')
        self.coarse = self.mesh.coarse
        self.mesh.permeability[:] = np.array([1, 1, 1]) # Diagonal permeability
        self.boundary_condition_type = boundary_condition_type
        self.coarse_config = coarse_config() # Access IMPRESS' internal class
        self.get_mesh_informations()
        self.number_coarse_volumes = len(self.coarse.elements) # Number of volumes from the coarse mesh
        self.number_volumes_local_problem = len(self.mesh.volumes)/(self.nx*self.ny*self.nz) # Number of fine scale volumes inside a coarse volume
        self.x = np.array([1,0,0])
        self.y = np.array([0,1,0])
        self.z = np.array([0,0,1])
        self.pressure_gradient = 1
        self.direction_string = np.array(['x', 'y', 'z'])
        self.directions_dictionary = {
            'x': self.x,
            'y': self.y,
            'z': self.z
            }

        # Set and solve local problems in x, y and z directions
        for i in self.direction_string:
            print('\nAssembly of local problems in {} direction...'.format(i))
            start = time.time()
            self.assembly_local_problem()
            end = time.time()
            print("\nThis step lasted {}".format(end-start))

            print('\nSetting boundary conditions in {} direction...'.format(i))
            start = time.time()
            self.set_boundary_conditions(i)
            end = time.time()
            print("\nThis step lasted {}".format(end-start))

            print('\nSolving local problems in {} direction...'.format(i))
            start = time.time()
            self.solve_local_problems()
            end = time.time()
            print("\nThis step lasted {}".format(end-start))

    def assembly_local_problem(self):
        """
        Assembly of local problems to generate the transmissibility matrix
        """

        for i in range(self.number_coarse_volumes):
            print("Assembly of local problem {}".format(i))
            faces_local_ids = self.coarse.elements[i].faces.internal # Local IDs from the internal faces from a coarse volume
            self.equivalent_permeability(i, faces_local_ids)
            adjacent_volumes = self.coarse.elements[i].faces.bridge_adjacencies(faces_local_ids, 2, 3).flatten() # Local IDs from both the neighbors from each of the internal faces
            neighbors_centers = np.reshape(self.coarse.elements[i].volumes.center[adjacent_volumes], newshape = (len(faces_local_ids), 6))
            neighbors_centers[:, 3:6] = neighbors_centers[:, 3:6]*(-1)
            centers_distance = np.linalg.norm((neighbors_centers[:, 0:3] + neighbors_centers[:, 3:6]), axis = 1); # Calculates the distante between the face's neighbors centers
            self.transmissibility = lil_matrix((int(self.number_volumes_local_problem), int(self.number_volumes_local_problem)), dtype = float)
            face_normal = self.coarse.elements[i].faces.normal[faces_local_ids]
            adjacent_volumes = np.reshape(adjacent_volumes, newshape = (len(faces_local_ids), 2))

            for j in range(len(faces_local_ids)):
                self.id1 = int(adjacent_volumes[j,0]) # ID of the first neighbor from the face
                self.id2 = int(adjacent_volumes[j,1]) # ID of the second neighbor from the face
                self.transmissibility[self.id1,self.id2] += self.mesh.coarse.elements[i].equivalent_permeability[faces_local_ids[j]]/centers_distance[j]
                self.transmissibility[self.id2,self.id1] += self.mesh.coarse.elements[i].equivalent_permeability[faces_local_ids[j]]/centers_distance[j]

            lil_matrix.setdiag(self.transmissibility,(-1)*self.transmissibility.sum(axis = 1))
            self.coarse.elements[i].transmissibility = self.transmissibility # Store transmissibility in a IMPRESS' variable

    def equivalent_permeability(self, coarse_volume, faces_local_ids):
        """
        Calculates the equivalent permeability of each internal face
        """
        faces_normal = self.coarse.elements[coarse_volume].faces.normal[faces_local_ids]
        permeability_column = np.arange(3)

        for i, j in zip(self.direction_string, permeability_column):
            direction = self.directions_dictionary.get(i)
            cross_product = np.cross(direction, faces_normal)
            norm_cross_product = np.linalg.norm(cross_product, axis = 1)
            self.index_faces_direction = np.isin(norm_cross_product, 0)
            faces_direction = faces_local_ids[self.index_faces_direction]
            adjacent_volumes = self.coarse.elements[coarse_volume].faces.bridge_adjacencies(faces_direction, 2, 3)
            global_id_adjacent_volumes = self.coarse.elements[coarse_volume].volumes.father_id[adjacent_volumes]
            permeability_direction = self.mesh.permeability[global_id_adjacent_volumes]
            permeability_direction = permeability_direction[:,j]
            permeability_direction = np.reshape(permeability_direction, newshape = (len(adjacent_volumes), 2))
            p1 = permeability_direction[:,0]
            p2 = permeability_direction[:,1]
            self.permeability_multiplication = np.multiply(p1, p2)
            self.permeability_sum = p1 + p2
            equivalent_permeability = 2*self.permeability_multiplication/self.permeability_sum
            self.global_ids_faces_direction = self.coarse.elements[coarse_volume].faces.father_id[faces_direction]
            self.coarse.elements[coarse_volume].equivalent_permeability[faces_direction] = equivalent_permeability
