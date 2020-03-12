"""
Module of local upscaling technique in structured tridimensional meshes using multiprocessing package to achieve greater efficiency
"""
import time
import numpy as np
from upscaling_procedures.local.parallel_local_problems import ParallelLocalProblems
from upscaling_procedures.local.local_upscaling import LocalUpscaling

class ParallelLocalUpscaling(ParallelLocalProblems, LocalUpscaling):

    def __init__(self, mesh_file = None, dataset = None):

        initial_time = time.time()

        print('\n##### Parallel local upscaling class initialized #####')

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

        # Upscale in parallel
        self.distribute_data()



        final_time = time.time()
        print("\nThe upscaling lasted {0}s".format(final_time-initial_time))

        def redistribute_pressure(self):
            pass

        def upscale_permeability_parallel(self):
            pass

        def upscale_porosity_parallel(self):
            pass

        def
