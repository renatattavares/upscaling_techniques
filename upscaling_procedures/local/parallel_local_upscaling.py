"""
Module of local upscaling technique in structured tridimensional meshes using multiprocessing package to achieve greater efficiency
"""
import time
import numpy as np
import multiprocessing as mp
from imex_integration.read_dataset import read_dataset
from upscaling_procedures.local.local_upscaling import LocalUpscaling
from upscaling_procedures.local.parallel_local_problems import ParallelLocalProblems
from impress.preprocessor.meshHandle.configTools.configClass import coarseningInit as coarse_config

class ParallelLocalUpscaling(ParallelLocalProblems, LocalUpscaling, ParallelFeatures):

    def __init__(self, mesh_file = None, dataset = None):

        initial_time = time.time()

        print('\n##### Parallel local upscaling class initialized #####')

        print('\n##### Treatment of local problems #####')

        if mesh_file is None:
            print('\nMesh informations will be accessed from {} dataset'.format(dataset))
            self.mode = 'integrated'
            self.mesh_file = 'mesh/generated_mesh.h5m'
            self.porosity, self.permeability = read_dataset(dataset, self.mesh_file)

        else:
            print('\nMesh informations will be set automatically')
            self.mode = 'auto'
            self.mesh_file = mesh_file

        # Read boundary condition chosen
        self.boundary_condition_type = 1

        # Preprocessing mesh with IMPRESS
        self.preprocess_mesh()

        # Setting variables and informations
        self.set_simulation_variables()
        self.set_coordinate_system()
        self.check_parallel_direction()
        self.get_mesh_informations(coarse_config())
        self.center_distance_walls()

        # Upscale in parallel
        self.distribute_data()
        self.create_processes()

        final_time = time.time()
        print("\nThe upscaling lasted {0}s".format(final_time-initial_time))

    def upscale_permeability(self, coarse_volume, pressure):

        area = 1
        i = coarse_volume
        effective_permeability = []

        for j in self.direction_string:
            direction = j
            if direction == 'x':
                local_wall = self.get_wall(i, 0)
                pressures = pressure[0]
                center_distance_walls = self.center_distance_walls_x[i]
            elif direction == 'y':
                local_wall = self.get_wall(i, 1)
                pressures = pressure[1]
                center_distance_walls = self.center_distance_walls_y[i]
            elif direction == 'z':
                local_wall = self.get_wall(i, 2)
                pressures = pressure[2]
                center_distance_walls = self.center_distance_walls_z[i]

            global_wall = self.coarse.elements[i].volumes.global_id[local_wall]
            local_adj = self.identify_adjacent_volumes_to_wall(i, local_wall)
            global_adj = self.coarse.elements[i].volumes.global_id[local_adj]
            center_wall = self.coarse.elements[i].volumes.center[local_wall]
            center_adj = self.coarse.elements[i].volumes.center[local_adj]
            pressure_wall = pressures[local_wall]
            pressure_adj = pressures[local_adj]
            permeability_wall = self.get_absolute_permeabilities(direction, global_wall)
            permeability_adj = self.get_absolute_permeabilities(direction, global_adj)

            flow_rate = ((2*np.multiply(permeability_wall,permeability_adj)/(permeability_wall+permeability_adj))*(pressure_wall-pressure_adj)/np.linalg.norm(center_wall - center_adj, axis = 1)).sum()

            effective_permeability.append(center_distance_walls*flow_rate/(area*self.number_faces_coarse_face))

        return effective_permeability

    def upscale_porosity(self):
        pass
