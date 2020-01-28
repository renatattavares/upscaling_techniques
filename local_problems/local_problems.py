"""
Module for treatment of local problems to apply local upscaling technique in structured tridimensional meshes
"""
import time
import numpy as np
#import xlsxwriter
import yaml
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from impress.preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as impress
from impress.preprocessor.meshHandle.configTools.configClass import coarseningInit as coarse_config
from local_problems.boundary_conditions import BoundaryConditionsLocalProblems


class LocalProblems:

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
        self.mesh.permeability[:] = np.array([1, 1, 1]) # Diagonal permeability tensor -> [Kx, 0, 0] -> [Kx, Ky, Kz]
                # [0, Ky, 0]
                # [0, 0, Kz]
        self.coarse = self.mesh.coarse
        self.boundary_condition_type = boundary_condition_type
        self.coarse_config = coarse_config() # Access IMPRESS' internal class
        self.get_mesh_informations()
        self.number_coarse_volumes = len(self.coarse.elements) # Number of volumes from the coarse mesh
        self.number_volumes_local_problem = len(self.mesh.volumes)/(self.nx*self.ny*self.nz) # Number of fine scale volumes inside a coarse volume
        self.x = np.array([1,0,0])
        self.y = np.array([0,1,0])
        self.z = np.array([0,0,1])
        self.pressure_gradient = 500
        self.system = np.array([self.x, self.y, self.z])

        # Initializing boundary conditions class
        bc = BoundaryConditionsLocalProblems(self, boundary_condition_type)

        # Assembly of local problems
        print('\nAssembly of local problems in x direction...')
        start = time.time()
        self.assembly_local_problem()
        end = time.time()
        print("\nThis step lasted {0}s".format(end-start))

        # Upscaled permeability in x direction
        print("\nSetting boundary conditions...")
        start = time.time()
        bc.set_boundary_conditions('x')
        end = time.time()
        print("\nThis step lasted {0}s".format(end-start))

        # Solve local problems in x direction
        print('\nSolving local problems...')
        start = time.time()
        self.solve_local_problems()
        end = time.time()
        print("\nThis step lasted {0}s".format(end-start))

        # Assembly of local problems
        print('\nAssembly of local problems in y direction...')
        start = time.time()
        self.assembly_local_problem()
        end = time.time()
        print("\nThis step lasted {0}s".format(end-start))
        # Upscaled permeability in y direction
        print("\nSetting boundary conditions...")
        start = time.time()
        bc.set_boundary_conditions('y')
        end = time.time()
        print("\nThis step lasted {0}s".format(end-start))
        # Solve local problems in y direction
        print('\nSolving local problems...')
        start = time.time()
        self.solve_local_problems()
        end = time.time()
        print("\nThis step lasted {0}s".format(end-start))

        # Assembly of local problems
        print('\nAssembly of local problems in z direction...')
        start = time.time()
        self.assembly_local_problem()
        end = time.time()
        print("\nThis step lasted {0}s".format(end-start))
        # Upscaled permeability in z direction
        print("\nSetting boundary conditions...")
        start = time.time()
        bc.set_boundary_conditions('z')
        end = time.time()
        print("\nThis step lasted {0}s".format(end-start))
        # Solve local problems in z direction
        print('\nSolving local problems...')
        start = time.time()
        self.solve_local_problems()
        end = time.time()
        print("\nThis step lasted {0}s".format(end-start))

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

    def assembly_local_problem(self):
        """
        Assembly of local problems to generate the transmissibility matrix
        """

        for i in range(self.number_coarse_volumes):
            print("Assembly of local problem {0}".format(i))
            faces_local_ids = self.coarse.elements[i].faces.internal # Local IDs from the internal faces from a coarse volume
            equivalent_permeability = self.equivalent_permeability(i, faces_local_ids)
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
                self.transmissibility[self.id1,self.id2] += equivalent_permeability[j]/centers_distance[j]
                self.transmissibility[self.id2,self.id1] += equivalent_permeability[j]/centers_distance[j]

            lil_matrix.setdiag(self.transmissibility,(-1)*self.transmissibility.sum(axis = 1))
            self.coarse.elements[i].transmissibility = self.transmissibility # Store transmissibility in a IMPRESS' variable

    def equivalent_permeability(self, coarse_volume, faces_local_ids):
        """
        Calculates the equivalent permeability of each internal face
        """
        faces_normal = self.coarse.elements[coarse_volume].faces.normal[faces_local_ids]
        equivalent_permeability = np.arange(len(faces_local_ids))

        # Testing faces for x direction
        cross_product = np.cross(self.x, faces_normal)
        norm_cross_product = np.linalg.norm(cross_product, axis = 1)
        index_faces_x_direction = np.isin(norm_cross_product, 0)
        faces_x_direction = faces_local_ids[index_faces_x_direction]
        adjacent_volumes = self.coarse.elements[coarse_volume].faces.bridge_adjacencies(faces_x_direction, 2, 3)
        global_id_adjacent_volumes_x = self.coarse.elements[coarse_volume].volumes.father_id[adjacent_volumes]
        permeability_x_direction = self.mesh.permeability[global_id_adjacent_volumes_x]
        permeability_x_direction = permeability_x_direction[:,0]
        permeability_x_direction = np.reshape(permeability_x_direction, newshape = (len(adjacent_volumes), 2))
        p1 = permeability_x_direction[:,0]
        p2 = permeability_x_direction[:,1]
        permeability_multiplication = np.multiply(p1, p2)
        permeability_sum = p1 + p2
        self.equivalent_permeability_x = 2*permeability_multiplication/permeability_sum

        # Testing faces for y direction
        cross_product = np.cross(self.y, faces_normal)
        norm_cross_product = np.linalg.norm(cross_product, axis = 1)
        index_faces_y_direction = np.isin(norm_cross_product, 0)
        faces_y_direction = faces_local_ids[index_faces_y_direction]
        adjacent_volumes = self.coarse.elements[coarse_volume].faces.bridge_adjacencies(faces_y_direction, 2, 3)
        global_id_adjacent_volumes_y = self.coarse.elements[coarse_volume].volumes.father_id[adjacent_volumes]
        permeability_y_direction = self.mesh.permeability[global_id_adjacent_volumes_y]
        permeability_y_direction = permeability_y_direction[:,1]
        permeability_y_direction = np.reshape(permeability_y_direction, newshape = (len(adjacent_volumes), 2))
        p1 = permeability_y_direction[:,0]
        p2 = permeability_y_direction[:,1]
        permeability_multiplication = np.multiply(p1, p2)
        permeability_sum = p1 + p2
        self.equivalent_permeability_y = 2*permeability_multiplication/permeability_sum

        # Testing faces for z direction
        cross_product = np.cross(self.z, faces_normal)
        norm_cross_product = np.linalg.norm(cross_product, axis = 1)
        index_faces_z_direction = np.isin(norm_cross_product, 0)
        faces_z_direction = faces_local_ids[index_faces_z_direction]
        adjacent_volumes = self.coarse.elements[coarse_volume].faces.bridge_adjacencies(faces_z_direction, 2, 3)
        global_id_adjacent_volumes_z = self.coarse.elements[coarse_volume].volumes.father_id[adjacent_volumes]
        permeability_z_direction = self.mesh.permeability[global_id_adjacent_volumes_z]
        permeability_z_direction = permeability_z_direction[:,2]
        permeability_z_direction = np.reshape(permeability_z_direction, newshape = (len(adjacent_volumes), 2))
        p1 = permeability_z_direction[:,0]
        p2 = permeability_z_direction[:,1]
        permeability_multiplication = np.multiply(p1, p2)
        permeability_sum = p1 + p2
        self.equivalent_permeability_z = 2*permeability_multiplication/permeability_sum
        #
        equivalent_permeability[index_faces_x_direction] = self.equivalent_permeability_x
        equivalent_permeability[index_faces_y_direction] = self.equivalent_permeability_y
        equivalent_permeability[index_faces_z_direction] = self.equivalent_permeability_z

        return equivalent_permeability

    def solve_local_problems(self):

        for i in range(self.number_coarse_volumes):
            print("Solving local problem of coarse volume {0}".format(i))
            transmissibility = lil_matrix.tocsr(self.coarse.elements[i].transmissibility)
            source = lil_matrix.tocsr(self.coarse.elements[i].source)
            global_id_volumes = self.coarse.elements[i].volumes.father_id[:]

            if np.array_equal(self.direction, self.x) == True:
                self.mesh.pressure_x[global_id_volumes] = spsolve(transmissibility,source)
            elif np.array_equal(self.direction, self.y) == True:
                self.mesh.pressure_y[global_id_volumes] = spsolve(transmissibility,source)
            elif np.array_equal(self.direction, self.z) == True:
                self.mesh.pressure_z[global_id_volumes] = spsolve(transmissibility,source)
