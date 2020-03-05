"""
Module of local upscaling technique in structured tridimensional meshes
"""
import time
import numpy as np
from upscaling_procedures.local.local_problems import LocalProblems

class LocalUpscaling:

    def __init__(self, mesh_file = None, boundary_condition_type = None, dataset = None):
        initial_time = time.time()

        if mesh_file is None:
            lp = LocalProblems(mesh_file = None, boundary_condition_type = 1, dataset = dataset)
        else:
            lp = LocalProblems(mesh_file = mesh_file, boundary_condition_type = 1, dataset = None)

        self.lp = lp

        # self.mesh = lp.mesh
        # self.coarse = lp.coarse
        # self.number_coarse_volumes = lp.number_coarse_volumes
        # self.number_faces_coarse_face = lp.number_faces_coarse_face
        # self.direction_string = lp.direction_string
        # self.directions_dictionary = lp.directions_dictionary
        # self.px = lp.pressure_x
        # self.py = lp.pressure_y
        # self.pz = lp.pressure_z
        # self.number_volumes_local_problem = lp.number_volumes_local_problem
        # self.p = lp.porosity
        # self.perm = lp.permeability
        # self.upscale_permeability()

        final_time = time.time()
        print("\nThe upscaling lasted {0}s".format(final_time-initial_time))

    def upscale_permeability(self):
        """
        It calculates the effective permeability of a coarse volume
        """
        for i in range(self.number_coarse_volumes):
            self.coarse_volume = i
            for j in np.array(['x']):
                self.direction = j
                wall = self.get_a_wall()
                self.wall = np.array(wall, dtype = int)
                self.adjacent_volumes = np.concatenate(self.coarse.elements[i].volumes.bridge_adjacencies(self.wall, 2, 3))
                repeated_volumes = np.isin(self.adjacent_volumes, self.wall, invert = True)
                self.adjacent_volumes = np.unique(self.adjacent_volumes[repeated_volumes])
                global_wall = self.coarse.elements[i].volumes.global_id[self.wall]
                global_adj = self.coarse.elements[i].volumes.global_id[self.adjacent_volumes]
                center_wall = self.coarse.elements[i].volumes.center[self.wall]
                center_adj = self.coarse.elements[i].volumes.center[self.adjacent_volumes]
                area = 1
                self.perm[:,0] = self.p + 5
                self.perm[:,1] = self.p + 5
                self.perm[:,2] = self.p + 5
                self.mesh.permeability[:] = self.perm

                if j == 'x':
                    pressure = self.px[int(i*self.number_volumes_local_problem):int((i*self.number_volumes_local_problem)+self.number_volumes_local_problem)]
                    permeability_wall = self.mesh.permeability[global_wall][:,0]
                    permeability_adj = self.mesh.permeability[global_adj][:,0]
                elif j == 'y':
                    pressure = self.py[int(i*self.number_volumes_local_problem):int((i*self.number_volumes_local_problem)+self.number_volumes_local_problem)]
                    permeability_wall = self.mesh.permeability[global_wall][:,1]
                    permeability_adj = self.mesh.permeability[global_adj][:,1]
                elif j == 'z':
                    pressure = self.pz[int(i*self.number_volumes_local_problem):int((i*self.number_volumes_local_problem)+self.number_volumes_local_problem)]
                    permeability_wall = self.mesh.permeability[global_wall][:,2]
                    permeability_adj = self.mesh.permeability[global_adj][:,2]

                pressure_wall = pressure[self.wall]
                pressure_adj = pressure[self.adjacent_volumes]
                self.flow_rate = (2*np.multiply(permeability_wall,permeability_adj)/(permeability_adj + permeability_wall)*(pressure_wall - pressure_adj)/np.linalg.norm(center_wall - center_adj, axis = 1)).sum()

                kef = 4*self.flow_rate/(area*self.number_faces_coarse_face)
                print(kef -1.7)
                global_id = self.coarse.elements[i].volumes.global_id[:]
                self.mesh.kefx[global_id] = kef - 1.7


    def upscale_porosity(self):
        """
        It calculates the effective porosity of a coarse volume
        """

        for i in range(self.number_coarse_volumes):
            global_id_volumes = self.coarse.elements[i].volumes.global_id[:]
            porosity = self.mesh.porosity[global_id_volumes]
            volumes = np.ones(len(porosity))
            effective_porosity = (np.multiply(porosity*volumes).sum())/volumes.sum()

    def get_a_wall(self):

        wall = np.array([])
        i = self.coarse_volume
        self.direction = self.directions_dictionary.get(self.direction)
        correct_volumes_group_1 = np.array([])
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

        global_ids_faces = self.coarse.elements[i].faces.global_id[:]
        index_group_1 = np.isin(global_ids_faces, group_1)
        index_group_2 = np.isin(global_ids_faces, group_2)
        local_ids_group_1 = np.where(index_group_1 == True)[0]
        local_ids_group_2 = np.where(index_group_2 == True)[0]

        correct_volumes_group_1 = self.coarse.elements[i].faces.bridge_adjacencies(local_ids_group_1, 2, 3).flatten()

        if np.array_equal(self.direction, np.array([1,0,0])) is True:
            wall = np.append(wall, correct_volumes_group_1)
        elif np.array_equal(self.direction, np.array([0,1,0])) is True:
            wall = np.append(wall, correct_volumes_group_1)
        elif np.array_equal(self.direction, np.array([0,0,1])) is True:
            wall = np.append(wall, correct_volumes_group_1)

        return wall
