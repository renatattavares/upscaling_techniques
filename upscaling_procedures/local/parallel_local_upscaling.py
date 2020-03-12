"""
Module of local upscaling technique in structured tridimensional meshes using multiprocessing package to achieve greater efficiency
"""
import time
import numpy as np
from upscaling_procedures.local.local_problems import LocalProblems
from upscaling_procedures.local.local_upscaling import LocalUpscaling

class ParallelLocalUpscaling(ParallelLocalProblems, LocalUpscaling):

    def __init__(self, mesh_file = None, dataset = None):


        if mesh_file is None:
            lp = LocalProblems(mesh_file = None, dataset = dataset)
        else:
            lp = LocalProblems(mesh_file = mesh_file, dataset = None)

        self.lp = lp

        final_time = time.time()
        print("\nThe upscaling lasted {0}s".format(final_time-initial_time))
