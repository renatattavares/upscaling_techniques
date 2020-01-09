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
        print('\nPre-processing mesh with IMPRESS...')
        start = time.time()
        self.mesh = preprocessor(mesh_file, dim = 3) # IMPRESS' object
        end = time.time()
        print("\nThis step lasted {0}s".format(end-start))

        # Setting variables
        print('\nAccessing coarsening informations from IMPRESS and setting important variables...')
        self.permeability = 1 # Inserir como propriedade do IMPRESS
        self.coarse = self.mesh.coarse
        self.boundary_condition_type = boundary_condition_type
        self.coarse_config = coarse_config() # Access IMPRESS' internal class
        self.get_coarse_informations()
        self.number_coarse_volumes = len(self.coarse.elements) # Number of volumes from the coarse mesh
        self.number_volumes_local_problem = len(self.mesh.volumes)/(self.nx*self.ny*self.nz) # Number of fine scale volumes inside a coarse volume
        self.x = np.array([1,0,0])
        self.y = np.array([0,1,0])
        self.z = np.array([0,0,1])

        # Assembly of local problems
        print('\nAssembly of local problems in x direction...')
        start = time.time()
        self.transmissibility = self.assembly_local_problem()
        end = time.time()
        print("\nThis step lasted {0}s".format(end-start))

        # Upscaled permeability in x direction
        print("\nSetting boundary conditions...")
        start = time.time()
        self.set_boundary_conditions('x')
        end = time.time()
        print("\nThis step lasted {0}s".format(end-start))

        # Solve local problems
        print('\nSolving local problems....\n')
        start = time.time()
        self.solve_local_problems()
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

    def assembly_local_problem(self):
        """
        Assembly of local problems to generate the transmissibility matrix
        """

        for i in range(self.number_coarse_volumes):
            print("Assembly of local problem {0}".format(i))
            local_ids = self.coarse.elements[i].faces.internal # Local IDs from the internal faces from a coarse volume
            global_ids = self.coarse.elements[i].faces.father_id[local_ids] # Global IDs from the internal faces from a coarse volume
            face_neighbors = self.coarse.elements[i].faces.bridge_adjacencies(local_ids, 2, 3) # Local IDs from both the neighbors from each of the internal faces
            first_neighbors_centers = self.coarse.elements[i].volumes.center[face_neighbors[:,0]] # Consult centers from the first column of the face's neighbors
            second_neighbors_centers = self.coarse.elements[i].volumes.center[face_neighbors[:,1]] # Consult centers from the second column of the face's neighbors
            centers_distance = np.linalg.norm((first_neighbors_centers - second_neighbors_centers), axis = 1) #Calculates the distante between the face's neighbors centers
            self.transmissibility = lil_matrix((int(self.number_volumes_local_problem), int(self.number_volumes_local_problem)), dtype = float)
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
        if self.direction.all() == self.x.all():
            self.perpendicular_direction_1 = self.y
            self.perpendicular_direction_2 = self.z

        elif self.direction.all() == self.y.all():
            self.perpendicular_direction_1 = self.x
            self.perpendicular_direction_2 = self.z

        elif self.direction.all() == self.z.all():
            self.perpendicular_direction_1 = self.x
            self.perpendicular_direction_2 = self.y

        boundary_conditions_dictionary.get(self.boundary_condition_type, "\nprint('Invalid boundary condition')")() # Execute the correct boundary condition function

    def fixed_constant_pressure(self):
        """
        Function to apply fixed constant pressure boundary condition, it returns a transmissibility and a source/sink matrix modified.
        """
        correct_volumes_group_1, correct_volumes_group_2 = self.identify_top_bottom_volumes()

        for i in range(self.number_coarse_volumes):
            volumes_group_1 = correct_volumes_group_1[i]
            volumes_group_2 = correct_volumes_group_2[i]
            self.coarse.elements[i].transmissibility[volumes_group_1] = 0
            self.coarse.elements[i].transmissibility[volumes_group_2] = 0
            self.coarse.elements[i].transmissibility[volumes_group_1, volumes_group_1] = 1
            self.coarse.elements[i].transmissibility[volumes_group_2, volumes_group_2] = 1
            self.coarse.elements[i].source = lil_matrix((int(self.number_volumes_local_problem), 1), dtype = 'float')
            self.coarse.elements[i].source[volumes_group_1] = 500
            self.coarse.elements[i].source[volumes_group_2] = 0

        print('\nFixed constant pressure boundary condition applied')

    def fixed_linear_pressure(self):
        """
        Function to apply fixed linear pressure boundary condition, it returns a transmissibility and a source/sink matrix modified.
        """

        print('\nFixed linear pressure boundary condition applied')

    def periodic_pressure(self):
        """
        Function to apply periodic pressure boundary condition, it returns a transmissibility and a source/sink matrix modified.
        """

        print('\nPeriodic pressure boundary condition applied')

    def identify_top_bottom_volumes(self):
        ########################## BAD DEFINITION ##########################
        correct_volumes_group_1 = np.zeros((self.number_coarse_volumes,25), dtype = int)
        correct_volumes_group_2 = np.zeros((self.number_coarse_volumes,25), dtype = int)
        ####################################################################

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

            adjacent_volumes_group_1 = self.coarse.elements[i].faces.bridge_adjacencies(local_ids_group_1, 2, 3)
            adjacent_volumes_group_1 = np.reshape(adjacent_volumes_group_1, newshape = (1,25))
            correct_volumes_group_1[i] = adjacent_volumes_group_1
            adjacent_volumes_group_2 = self.coarse.elements[i].faces.bridge_adjacencies(local_ids_group_2, 2, 3)
            adjacent_volumes_group_2 = np.reshape(adjacent_volumes_group_2, newshape = (1,25))
            correct_volumes_group_2[i] = adjacent_volumes_group_2

        return correct_volumes_group_1, correct_volumes_group_2

    def identify_side_volumes(self):
        ########################## BAD DEFINITION ##########################
        correct_volumes_group_1 = np.zeros((self.number_coarse_volumes,25), dtype = int)
        correct_volumes_group_2 = np.zeros((self.number_coarse_volumes,25), dtype = int)
        correct_volumes_group_3 = np.zeros((self.number_coarse_volumes,25), dtype = int)
        correct_volumes_group_4 = np.zeros((self.number_coarse_volumes,25), dtype = int)
        ####################################################################

        for i in range(self.number_coarse_volumes):
            boundary_faces = self.coarse.elements[i].faces.boundary # Local IDs of boundary faces of a coarse volume
            normal_boundary_faces = self.coarse.elements[i].faces.normal[boundary_faces] # Normal vector of boundary faces of a coarse volumes
            direction_vector = np.ndarray(shape = np.shape(normal_boundary_faces), dtype = float)
            direction_vector = np.full_like(direction_vector, self.perpendicular_direction_1) # Vectorization

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

            global_ids_faces = self.coarse.elements[i].faces.father_id[:]
            index_group_1 = np.isin(global_ids_faces, group_1)
            index_group_2 = np.isin(global_ids_faces, group_2)
            local_ids_group_1 = np.where(index_group_1 == True)
            local_ids_group_1 = np.where(index_group_2 == True)

            adjacent_volumes_group_1 = self.coarse.elements[i].faces.bridge_adjacencies(local_ids_group_1, 2, 3)
            correct_volumes_group_1[i] = adjacent_volumes_group_1
            adjacent_volumes_group_2 = self.coarse.elements[i].faces.bridge_adjacencies(local_ids_group_2, 2, 3)
            correct_volumes_group_2[i] = adjacent_volumes_group_2

    def solve_local_problems(self):

        for i in range(self.number_coarse_volumes):
            print("Solving local problem of coarse volume {0}".format(i))
            transmissibility = lil_matrix.tocsr(self.coarse.elements[i].transmissibility)
            source = lil_matrix.tocsr(self.coarse.elements[i].source)
            global_id_volumes = self.coarse.elements[i].volumes.father_id[:]
            self.mesh.pressure[global_id_volumes] = spsolve(transmissibility,source)
