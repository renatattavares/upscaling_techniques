"""
Module of local upscaling technique in structured tridimensional meshes
"""
import time
import numpy as np
from upscaling_procedures.local.local_problems import LocalProblems

class LocalUpscaling(LocalProblems):

    def __init__(self, mesh_file = None, dataset = None):
        initial_time = time.time()

        super().__init__(mesh_file, dataset)
        self.center_distance_walls()
        self.upscale_permeability()

        final_time = time.time()
        print("\nThe upscaling lasted {0}s".format(final_time-initial_time))


    def upscale_permeability(self):
        """
        It calculates the effective permeability of a coarse volume
        """
        area = 1
        for i in range(self.number_coarse_volumes):
            for j in self.direction_string:
                direction = j

                if direction == 'x':
                    wall = self.coarse_face_x[i]
                    pressures = self.pressure_x
                    center_distance_walls = self.center_distance_walls_x[i]
                elif direction == 'y':
                    wall = self.coarse_face_y[i]
                    pressures = self.pressure_y
                    center_distance_walls = self.center_distance_walls_y[i]
                elif direction == 'z':
                    wall = self.coarse_face_z[i]
                    pressures = self.pressure_z
                    center_distance_walls = self.center_distance_walls_z[i]

                global_wall = self.coarse.elements[i].volumes.global_id[wall]
                local_adj = self.identify_adjacent_volumes_to_wall(i, wall)
                global_adj = self.coarse.elements[i].volumes.global_id[local_adj]
                center_wall = self.coarse.elements[i].volumes.center[wall]
                center_adj = self.coarse.elements[i].volumes.center[local_adj]
                self.pressure_wall = pressures[global_wall]
                self.pressure_adj = pressures[global_adj]
                self.permeability_wall = self.get_absolute_permeabilities(direction, global_wall)
                self.permeability_adj = self.get_absolute_permeabilities(direction, global_adj)

                # flow_rate = (2*np.multiply(permeability_wall,permeability_adj)/(permeability_wall+permeability_adj))*(pressure_wall-pressure_adj)/np.linalg.norm(center_wall - center_adj, axis = 1).sum()
                #
                # effective_permeability = center_distance_walls[i]*flow_rate/(area*self.number_faces_coarse_face)
                #
                # print(effective_permeability)

    def upscale_porosity(self):
        """
        It calculates the effective porosity of a coarse volume
        """

        for i in range(self.number_coarse_volumes):
            global_id_volumes = self.coarse.elements[i].volumes.global_id[:]
            porosity = self.mesh.porosity[global_id_volumes]
            volumes = np.ones(len(porosity))
            effective_porosity = (np.multiply(porosity*volumes).sum())/volumes.sum()

    def identify_adjacent_volumes_to_wall(self, coarse_volume, wall):

        adjacent_volumes = np.concatenate(self.coarse.elements[coarse_volume].volumes.bridge_adjacencies(wall, 2, 3))
        repeated_volumes = np.isin(adjacent_volumes, wall, invert = True)
        adjacent_volumes = np.unique(adjacent_volumes[repeated_volumes].flatten())

        return adjacent_volumes

    def get_absolute_permeabilities(self, direction, global_ids):

        if direction == 'x':
            permeability = self.mesh.permeability[global_ids][0]
        elif direction == 'y':
            permeability = self.mesh.permeability[global_ids][1]
        elif direction == 'z':
            permeability = self.mesh.permeability[global_ids][2]

        return permeability

    def center_distance_walls(self):

        self.center_distance_walls_x = np.array([])
        self.center_distance_walls_y = np.array([])
        self.center_distance_walls_z = np.array([])

        for i in range(self.number_coarse_volumes):
            min_coord = self.coarse.elements[i].volumes.center[:].min(axis = 0)
            max_coord = self.coarse.elements[i].volumes.center[:].min(axis = 0)
            self.center_distance_walls_x = np.append(self.center_distance_walls_x, (min_coord[0], max_coord[0]))
            self.center_distance_walls_y = np.append(self.center_distance_walls_x, (min_coord[1], max_coord[1]))
            self.center_distance_walls_z = np.append(self.center_distance_walls_x, (min_coord[2], max_coord[2]))
