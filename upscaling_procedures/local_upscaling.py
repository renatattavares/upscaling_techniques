"""
Module of local upscaling technique in structured tridimensional meshes
"""
import time
import numpy as np
import xlsxwriter
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

class LocalUpscaling:

    def __init__(self, preprocessor, coarse_config, mesh_file = None, boundary_condition_type = None):

        print('\n##### Local upscaling class initialized #####')

        # Preprocessing mesh with IMPRESS
        print('\nPre-processing mesh with IMPRESS')
        start = time.time()
        self.mesh = preprocessor(mesh_file, dim = 3) # IMPRESS' object
        end = time.time()
        print("\nThis step lasted {0}s".format(end-start))

        # Setting variables
        print('\nAccessing coarsening informations from IMPRESS...')
        self.permeability = 1 # Inserir como propriedade do IMPRESS
        self.coarse = self.mesh.coarse
        self.boundary_condition_type = boundary_condition_type
        self.coarse_config = coarse_config() # Access IMPRESS' internal class
        self.get_coarse_informations()
        self.number_coarse_volumes = len(self.coarse.elements) # Number of volumes from the coarse mesh
        self.number_volumes_local_problem = len(self.mesh.volumes)/(self.nx*self.ny*self.nz) # Number of fine scale volumes inside a coarse volume
        self.get_coarse_centers()

        # Coordinate System
        print("Setting coordinates system...")
        self.x = np.array([1,0,0])
        self.y = np.array([0,1,0])
        self.z = np.array([0,0,1])

        # Assembly of local problems
        print('\nAssembly of local problems in x direction')
        start = time.time()
        self.transmissibility = self.assembly_local_problem()
        end = time.time()
        print("\nThis step lasted {0}s".format(end-start))

        # Upscaled permeability in x direction
        print("\nSetting boundary conditions")
        start = time.time()
        self.get_coarse_neighbors()
        #self.set_boundary_conditions('x')
        end = time.time()
        print("\nThis step lasted {0}s".format(end-start))

        # Solve local problems
        print('\nSolving local problems\n')
        start = time.time()
        #self.solve_local_problems()
        end = time.time()
        print("\nThis step lasted {0}s".format(end-start))

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
            print("Local problem {0}".format(i))
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
            self.coarse.elements[i].transmissibility = self.transmissibility # Store transmissibility in a IMPRESS' variable

    def set_boundary_conditions(self, direction):
        """
        Indicates which function must be executed to set boundary conditions acording to the option informed by the user.
        """

        directions_dictionary = {
            'x': self.x,
            'y': self.y,
            'z': self.z
            }

        boundary_conditions_dictionary = {
            1: self.fixed_constant_pressure, # Fixed constant pressure
            2: self.fixed_linear_pressure,   # Fixed linear pressure
            3: self.periodic_pressure        # Periodic pressure
            }

        self.direction = directions_dictionary.get(direction) # Get the direction given

        boundary_conditions_dictionary.get(self.boundary_condition_type, "\nprint('Invalid boundary condition')")() # Execute the correct boundary condition function

    def fixed_constant_pressure(self):
        """
        Must return a mesh, a transmissibility and a source/sink matrix with boundary conditions set.
        """
        self.correct_volumes, self.correct_faces = self.get_correct_volumes_and_separate_faces()

        for i in range(self.number_coarse_volumes):
            self.correct_volumes_local_ids = self.correct_volumes_local_ids[i]

            self.coarse.elements[i].transmissibility[self.correct_volumes_local_ids] = 0
            #self.coarse.elements[i].transmissibility[self.correct_volumes_local_ids, self.correct_volumes_local_ids] = 1
            #self.coarse.elements[i].source = lil_matrix((int(self.number_volumes_local_problem), 1), dtype = 'float')
            #self.coarse.elements[i].source[self.correct_volumes_local_ids] = 500

        print('\nFixed constant pressure boundary condition applied')

    def get_correct_volumes_and_separate_faces(self):

        correct_volumes = np.zeros((self.number_coarse_volumes,50), dtype = int)

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
            adjacent_volumes = self.coarse.elements[i].faces.bridge_adjacencies(correct_faces, 2, 3)
            correct_volumes[i] = adjacent_volumes

            # Separate faces in two groups
            global_ids_correct_faces = self.coarse.elements[i].faces.global_ids[correct_faces]
            interface_coarse_face_id = self.coarse.iface_neighbors(i)[1]

            for j in range(interface_coarse_face_id):
                interface_faces = self.coarse.elements[i].interface_faces[j]
                verify = np.any(np.isin(interface_faces, global_ids_correct_faces[0]))
                if verify == True:
                    group_1 = np.isin(interface_faces, global_ids_correct_faces)
                    index = np.where(group_1 == True)
                    group_1 = interface_faces[index]
                    break

            group_2 = np.isin(global_ids_correct_faces, group_1)
            index = np.where(group_2 == True)
            group_2 = global_ids_correct_faces[index]

        print('\nCoarse neighbors defined')

        return correct_volumes

    def solve_local_problems(self):

        for i in range(self.number_coarse_volumes):
            print("Solving local problem of coarse volume {0}".format(i))
            transmissibility = lil_matrix.tocsr(self.coarse.elements[i].transmissibility)
            source = lil_matrix.tocsr(self.coarse.elements[i].source)
            self.coarse.elements[i].coarse_pressure = spsolve(transmissibility,source)

    def upscaled_permeability(self):
        pass

    def upscaled_transmissibility(self):
        pass

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
