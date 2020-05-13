import time
import numpy as np
import multiprocessing as mp
from imex_integration.read_dataset import read_dataset
from upscaling_procedures.local.local_upscaling import LocalUpscaling
from upscaling_procedures.local.parallel_local_problems import ParallelLocalProblems
from impress.preprocessor.meshHandle.configTools.configClass import coarseningInit as coarse_config

class ParallelFeatures:

    def get_wall(self, coarse_volume, direction):

        i = coarse_volume
        wall = np.zeros((1, self.number_faces_coarse_face), dtype = int)

        boundary_faces = self.coarse.elements[i].faces.boundary # Local IDs of boundary faces of a coarse volume
        global_ids_faces = self.coarse.elements[i].faces.global_id[boundary_faces]
        parallel_direction = self.mesh.parallel_direction[global_ids_faces]
        index_faces_direction = np.isin(parallel_direction, direction)
        index_faces_direction = np.where(index_faces_direction == True)[0]
        correct_faces = boundary_faces[index_faces_direction]

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

        global_ids_faces = self.coarse.elements[i].faces.global_id[:]
        index_group_1 = np.isin(global_ids_faces, group_1)
        local_ids_group_1 = np.where(index_group_1 == True)[0]

        wall = self.coarse.elements[i].faces.bridge_adjacencies(local_ids_group_1, 2, 3).flatten()

        return wall

    def upscale_permeability_parallel(self, coarse_volumes, queue):

        ep = []

        for cv in coarse_volumes:
            p = self.solve_local_problems(cv)
            ep.append(self.upscale_permeability(cv, p))
            print('Done with coarse volume {}'.format(cv))

        queue.put(ep)

    def create_processes(self):

        effective_permeabilities = []
        process_number = len(self.distribution)

        queues = [mp.Queue() for i in range(process_number)]
        processes = [mp.Process(target = self.upscale_permeability_parallel, args = (i,q,)) for i, q in zip(self.distribution, queues)]

        for p in processes:
            p.start()

        for p, q in zip(processes, queues):
            effective_permeabilities.append(q.get())
            p.join()

        self.effective_permeability = effective_permeabilities

    def print_results(self):

        temp_dist = np.array([], dtype = int)
        temp_perm = np.array([])

        for i in self.distribution:
            temp = np.array([(result) for result in i])
            temp_dist = np.append(temp_dist, temp)

        for i in self.effective_permeability:
            temp = np.array([(result[0]) for result in i])
            temp_perm = np.append(temp_perm, temp)

        for problem, perm in zip(temp_dist, temp_perm):
            self.coarse.elements[problem].kefx[:] = perm
