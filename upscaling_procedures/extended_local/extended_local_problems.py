"""
Module for treatment of extended local problems to apply local upscaling technique in structured tridimensional meshes
"""
import time
import numpy as np
from upscaling_procedures.local.local_problems import LocalProblems
from upscaling_procedures.extended_local.extended_assembly import ExtendedAssembly
from impress.preprocessor.meshHandle.configTools.configClass import coarseningInit as coarse_config

class ExtendedLocalProblems(LocalProblems, ExtendedAssembly):

    def __init__(self, mesh_file = None, dataset = None):
        initial_time = time.time()

        print('\n##### ExtendedLocalProblems class initialized #####')

        if mesh_file is None:
            self.mesh_file = 'mesh/dataset_mesh.h5m'
            porosity, permeability, self.number_elements, self.length_elements = read_dataset(dataset, self.mesh_file)
            # Preprocessing mesh with IMPRESS
            self.preprocess_mesh()
            self.mesh.porosity[:] = porosity
            self.mesh.permeability[:] = permeability

        else:
            self.mode = 'auto'
            self.mesh_file = mesh_file
            # Preprocessing mesh with IMPRESS
            self.preprocess_mesh()
            self.set_simulation_variables()

        # Setting variables and informations
        self.boundary_condition_type = 1 # Fixed pressure
        self.set_coordinate_system()

        # Preparing to solve local problems
        self.check_parallel_direction()
        self.get_mesh_informations(coarse_config()) # Remenber number_volumes_local_problem must be corrected (calculate it for each problem)
        self.transmissibility = self.assembly_extended_local_problem(0)

    def volumes_extended_local_problem(self, coarse_volume):
        """
        Define temporary local ids for an extended local problem
        """
        coarse_neighbors = self.coarse.inode_neighbors(coarse_volume)[0]
        delete_default_value, index, discard = np.intersect1d(coarse_neighbors, self.number_coarse_volumes, return_indices = True)

        if delete_default_value.size != 0:
            coarse_neighbors = np.delete(coarse_neighbors, index)

        global_ids_extended_local_problem = []
        global_ids_extended_local_problem.append(self.coarse.elements[coarse_volume].volumes.father_id[:])

        for neighbor in coarse_neighbors:
            global_ids_extended_local_problem.append(self.coarse.elements[neighbor].volumes.father_id[:])

        global_ids_extended_local_problem = np.array(global_ids_extended_local_problem).flatten()
        local_ids_extended_local_problem = np.arange(len(global_ids_extended_local_problem), dtype = int)

        return global_ids_extended_local_problem, local_ids_extended_local_problem

    def faces_extended_local_problem(self, coarse_volume):

        volumes_global_ids, volumes_local_ids = self.volumes_extended_local_problem(coarse_volume)
        faces_extended_local_problem = np.unique(self.mesh.volumes.adjacencies[volumes_global_ids].flatten())
        faces_centers = self.mesh.faces.center[faces_extended_local_problem]
        boundary_indexes = []

        for direction in self.direction_string:
            direction_number = self.directions_numbers.get(direction)
            coords = faces_centers[:, direction_number]
            greatest = np.amax(coords)
            smallest = np.amin(coords)
            boundary_indexes.append(np.where(coords == greatest)[0])
            boundary_indexes.append(np.where(coords == smallest)[0])

        boundary_indexes = np.array(boundary_indexes).flatten()
        boundary_faces_extended_local_problem = faces_extended_local_problem[boundary_indexes]
        internal_faces_extended_local_problem = np.delete(faces_extended_local_problem, boundary_indexes)

        return faces_extended_local_problem, boundary_faces_extended_local_problem, internal_faces_extended_local_problem
