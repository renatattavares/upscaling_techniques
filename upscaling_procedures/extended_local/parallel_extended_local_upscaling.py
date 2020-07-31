"""
Module of local upscaling technique in structured tridimensional meshes using multiprocessing package to achieve greater efficiency
"""
import time
import yaml
import numpy as np
import multiprocessing as mp
from imex_integration.read_dataset import read_dataset
from upscaling_procedures.local.parallel_local_upscaling import ParallelLocalUpscaling
from upscaling_procedures.extended_local.extended_local_upscaling import ExtendedLocalUpscaling
from impress.preprocessor.meshHandle.configTools.configClass import coarseningInit as coarse_config

class ParallelExtendedLocalUpscaling(ExtendedLocalUpscaling):

    def __init__(self, mesh_file = None, dataset = None):

        super().__init__(mesh_file, dataset)

        # # Check is refinement is required
        # refine = self.check_if_refinement_is_required()
        #
        # if refine is True:
        #     self.dont_upscale = self.identify_coarse_volumes_in_refinement_regions() # These coarse volumes shouldn't be upscaled
        # else:
        #     self.dont_upscale = np.array([])

        self.dont_upscale = np.array([])

        # Upscale in parallel
        self.distribute_data()
        self.create_processes()

    def upscale_permeability_porosity(self, cv):

        print('Upscaling of local problem {}'.format(cv))
        volume = self.length_elements[0]*self.length_elements[1]*self.length_elements[2]
        general_transmissibility = self.assembly_extended_local_problem(cv)
        info = []
        volumes_global_ids, volumes_local_ids = self.volumes_extended_local_problem(cv)
        total_volume = volume*len(volumes_global_ids)
        porosity = self.mesh.porosity[volumes_global_ids]
        effective_porosity = (volume*porosity).sum()/total_volume

        for direction in self.direction_string:
            #print('in {} direction'.format(direction))
            self.direction = direction
            transmissibility, source, local_wall = self.set_boundary_conditions(general_transmissibility, cv)
            pressures = self.solver(transmissibility, source)

            direction_number = self.directions_numbers.get(direction)
            center_distance_walls = self.center_distance_walls[cv][direction_number]
            area = self.areas[direction_number]
            global_wall = volumes_global_ids[local_wall]
            global_adj = self.identify_adjacent_volumes_to_wall(cv, global_wall)
            center_wall = self.mesh.volumes.center[global_wall]
            center_adj = self.mesh.volumes.center[global_adj]
            local_adj = self.global_to_local_id(volumes_global_ids, global_adj)
            pressure_wall = pressures[local_wall]
            pressure_adj = pressures[local_adj]
            pw = self.mesh.permeability[global_wall]
            permeability_wall = pw[:, direction_number]
            pa = self.mesh.permeability[global_adj]
            permeability_adj = pa[:, direction_number]
            flow_rate = ((((2*np.multiply(permeability_wall,permeability_adj)/(permeability_wall+permeability_adj))*(pressure_wall-pressure_adj))*self.areas[direction_number])/np.linalg.norm(center_wall - center_adj, axis = 1)).sum()
            info.append(center_distance_walls*flow_rate/(area*len(global_wall)))

        info.append(effective_porosity)

        return info

    def identify_adjacent_volumes_to_wall(self, coarse_volume, global_wall):

        adjacent_volumes = np.concatenate(self.mesh.volumes.bridge_adjacencies(global_wall, 2, 3))
        volumes_global_ids, volumes_local_ids = self.volumes_extended_local_problem(coarse_volume)
        volumes_inside_extended_problem = np.isin(adjacent_volumes, volumes_global_ids)
        volumes_inside_extended_problem = adjacent_volumes[volumes_inside_extended_problem]
        repeated_volumes = np.isin(volumes_inside_extended_problem, global_wall, invert = True)
        adjacent_volumes = np.unique(volumes_inside_extended_problem[repeated_volumes].flatten())

        return adjacent_volumes

    def upscale_permeability_porosity_parallel(self, coarse_volumes, queue):

        info = []

        for cv in coarse_volumes:
            info.append(self.upscale_permeability_porosity(cv))

        queue.put(info)

    def distribute_data(self):

        core_number = mp.cpu_count() # Cores in the computer
        number_problems_in_process = int(self.number_coarse_volumes/core_number) # Average number of local problems to be solved by one process
        rest = self.number_coarse_volumes % core_number
        distribution = []
        coarse_volumes = np.arange(self.number_coarse_volumes)

        if self.dont_upscale.size is not 0:
            inter, comm1, comm2 = np.intersect1d(coarse_volumes, self.dont_upscale, return_indices = True)
            coarse_volumes = np.delete(coarse_volumes, comm1)
        else:
            pass

        for i in range(core_number):
            if i is core_number-1:
                distribution.append(coarse_volumes[number_problems_in_process*i:number_problems_in_process*(i+1)+rest])
            else:
                distribution.append(coarse_volumes[number_problems_in_process*i:number_problems_in_process*(i+1)])

        self.distribution = distribution

    def create_processes(self):

        self.info = []
        process_number = len(self.distribution)

        queues = [mp.Queue() for i in range(process_number)]
        processes = [mp.Process(target = self.upscale_permeability_porosity_parallel, args = (i,q,)) for i, q in zip(self.distribution, queues)]

        for p in processes:
            p.start()

        for p, q in zip(processes, queues):
            self.info.append(q.get())
            p.join()

        self.final_time = time.time()
