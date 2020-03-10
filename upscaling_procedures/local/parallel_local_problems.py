"""
Module for treatment of local problems to apply local upscaling technique in structured tridimensional meshes using multiprocessing package to achieve greater efficiency
"""
import time
import numpy as np
from multiprocessing import Process
from upscaling_procedures.local.local_problems import LocalProblems
from impress.preprocessor.meshHandle.configTools.configClass import coarseningInit as coarse_config

class ParallelLocalProblems(LocalProblems):

    def __init__(self, mesh_file = None, dataset = None):

        print('\n##### Treatment of local problems #####')

        if mesh_file is None:
            print('\nMesh informations will be accessed from {} dataset'.format(dataset))
            self.mode = 'integrated'
            self.porosity, self.permeability = read_dataset(dataset)
            self.mesh_file = 'generated_mesh.h5m'

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

        self.create_processes()

        # Solve local problems
        #self.solve_local_problems()

    def distribute_data(self):
        pass


    def create_processes(self):
        p = Process(target = self.set_coordinate_system, args = (self,))
