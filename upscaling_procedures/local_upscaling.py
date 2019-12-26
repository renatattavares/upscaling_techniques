"""
Module of local upscaling technique in structured tridimensional meshes
"""
import time
import numpy as np
import xlsxwriter
from scipy.sparse import lil_matrix
#from scipy.sparse.linalg import spsolve

class LocalUpscaling:

    def __init__(self, preprocessor, coarse_config, mesh_file = None, boundary_condition_type = None):

        print('\n##### Local upscaling class initialized #####')

        # Preprocessing mesh with IMPRESS
        print('\nPre-processing mesh with IMPRESS')
        start = time.time()
        self.mesh = preprocessor(mesh_file, dim = 3) # IMPRESS' object
        end = time.time()
        print("\nThe pre-processing step lasted {0}s".format(end-start))

        # Setting variables
        print('\nAccessing coarsening informations from IMPRESS')
        self.permeability = 1 # Inserir como propriedade do IMPRESS
        self.coarse = self.mesh.coarse
        self.boundary_condition_type = boundary_condition_type
        self.coarse_config = coarse_config() # Access IMPRESS' internal class
        self.get_coarse_informations()
        self.number_coarse_volumes = len(self.coarse.elements) # Number of volumes from the coarse mesh
        self.number_volumes_local_problem = len(self.mesh.volumes)/(self.nx*self.ny*self.nz) # Number of fine scale volumes inside a coarse volume
        self.get_coarse_centers()

        # Coordinate System
        print('\nSetting coordinates system')
        self.x = np.array([1,0,0])
        self.y = np.array([0,1,0])
        self.z = np.array([0,0,1])

        # Assembly of local problems
        print('\nAssembly of local problems in x direction')
        start = time.time()
        self.transmissibility = self.assembly_local_problem()
        end = time.time()
        print("\nThe assembly lasted {0}s".format(end-start))

        # Upscaled permeability in x direction
        # self.transmissibility, self.flux = self.set_boundary_conditions('x')

    def get_coarse_informations(self):
        """
        Access coarsening informations given to IMPRESS.
        """
        self.tree = self.coarse_config.tree['Simple']
        self.nx = self.tree['nx']
        self.ny = self.tree['ny']
        self.nz = self.tree['nz']

        print('\nCoarsening informations accessed')

    def get_coarse_centers(self):
        """
        Calculates the center of a coarse volume and stores it in a IMPRESS variable called coarse_center
        """
        for i in range(self.number_coarse_volumes):
            coords = self.coarse.elements[i].nodes.coords[:]
            x_coords = coords[:,0]
            y_coords = coords[:,1]
            z_coords = coords[:,2]

            xmax = x_coords.max()
            ymax = y_coords.max()
            zmax = z_coords.max()

            xmin = x_coords.min()
            ymin = y_coords.min()
            zmin = z_coords.min()

            x_center = xmin + (xmax-xmin)/2
            y_center = ymin + (ymax-ymin)/2
            z_center = zmin + (zmax-zmin)/2

            self.coarse.elements[i].coarse_center[:] = np.array([x_center, y_center, z_center])

        print('\nCoarse volume centers calculated')

    def assembly_local_problem(self):
        """
        Assembly of local problems to generate the transmissibility matrix
        """

        for i in range(self.number_coarse_volumes):
            local_ids = self.coarse.elements[i].faces.internal # Local IDs from the internal faces from a coarse volume
            global_ids = self.coarse.elements[i].faces.father_id[local_ids] # Global IDs from the internal faces from a coarse volume
            face_neighbors = self.coarse.elements[i].faces.bridge_adjacencies(local_ids, 2, 3) # Local IDs from both the neighbors from each of the internal faces
            first_neighbors_centers = self.coarse.elements[i].volumes.center[face_neighbors[:,0]] # Consult centers from the first column of the face's neighbors
            second_neighbors_centers = self.coarse.elements[i].volumes.center[face_neighbors[:,1]] # Consult centers from the second column of the face's neighbors
            centers_distance = np.linalg.norm((first_neighbors_centers - second_neighbors_centers), axis = 1) #Calculates the distante between the face's neighbors centers
            self.transmissibility = lil_matrix((int(self.number_volumes_local_problem), int(self.number_volumes_local_problem)), dtype = 'float')
            face_normal = self.coarse.elements[i].faces.normal[local_ids]

            for j in range(len(local_ids)):
                self.id1 = int(face_neighbors[j,0]) # ID of the first neighbor from the face
                self.id2 = int(face_neighbors[j,1]) # ID of the second neighbor from the face
                self.transmissibility[self.id1,self.id2] += self.permeability/centers_distance[j]
                self.transmissibility[self.id2,self.id1] += self.permeability/centers_distance[j]

            lil_matrix.setdiag(self.transmissibility,(-1)*self.transmissibility.sum(axis = 1))

        return self.transmissibility


    def get_coarse_neighbors(self):

        self.correct_neighbors = np.zeros((self.number_coarse_volumes, 2))

        for i in range(self.number_coarse_volumes):
            neighbors = self.coarse.elements[i].faces.coarse_neighbors
            boundary_element = np.isin(neighbors, self.number_coarse_volumes)
            index = np.where(boundary_element == True)
            neighbors = np.delete(neighbors, index)

            for j in neighbors:
                center = self.coarse.elements[j].coarse_center[0]
                vector = self.coarse.elements[i].coarse_center[0] - center
                pv = np.cross(vector, self.direction)

                if np.linalg.norm(pv) == 0:
                    if self.correct_neighbors[i,0] == 0:
                        self.correct_neighbors[i,0] = j
                    else:
                        self.correct_neighbors[i,1] = j

        print('\nCoarse neighbors defined')


    def set_boundary_conditions(self, direction):
        """
        Indicates which function must be executed to set boundary conditions acording to the option informed by the user.
        """

        self.direction = direction

        self.directions = {
            'x': self.x,
            'y': self.y,
            'z': self.z
            }

        self.boundary_conditions = {
            1: 'self.fixed_constant_pressure()', # Fixed constant pressure
            2: 'self.fixed_linear_pressure()',   # Fixed linear pressure
            3: 'self.periodic_pressure()'        # Periodic pressure
            }

        self.direction = self.directions.get(self.direction)

        self.bc = self.boundary_conditions.get(self.boundary_condition_type, "\nprint('Invalid boundary condition')")

        self.get_coarse_neighbors()

        exec(self.bc)


    def fixed_constant_pressure(self):
        """
        Must return a mesh with boundary conditions set in the specific given volume
        """
        self.neighbors = np.amax(self.correct_neighbors, axis = 1)

        for i in range(self.number_coarse_volumes):
            coarse_neighbors = self.coarse.iface_neighbors(i)[0]
            coarse_faces = self.coarse.iface_neighbors(i)[1]
            boundary_element = np.isin(coarse_neighbors, self.neighbors[i])
            index = np.where(boundary_element == True)[0]
            interface = coarse_faces[index][0]
            faces = self.coarse.interfaces_faces[interface.item()]
            adj_volumes = self.mesh.faces.bridge_adjacencies(faces,2,3)
            adj_volumes = np.concatenate(adj_volumes, axis = 0)
            adj_volumes = np.unique(adj_volumes)
            fine_volumes = self.coarse.elements[i].volumes.father_id[:]
            correct_volumes = np.isin(fine_volumes, adj_volumes)
            index =  np.where(correct_volumes == True)
            correct_volumes = fine_volumes[index]
            self.mesh.pressure[correct_volumes] = 1000

        print('Fixed constant pressure boundary condition applied')


    def fixed_linear_pressure(self):
        """
        Must return a mesh with boundary conditions set in the specific given volume
        """

        print('\nFixed linear pressure boundary condition applied')

    def periodic_pressure(self):
        """
        Must return a mesh with boundary conditions set in the specific given volume
        """

        print('\nPeriodic pressure boundary condition applied')

    def solver(self):
        pass

    def upscaled_permeability(self):
        pass

    def upscaled_transmissibility(self):
        pass
