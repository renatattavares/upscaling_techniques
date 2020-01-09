"""
Module of local upscaling technique in structured tridimensional meshes
"""
from impress.preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as impress
from impress.preprocessor.meshHandle.configTools.configClass import coarseningInit as coarse_config
from local_problems import LocalProblems
import time
import numpy as np

class LocalUpscaling:

    def __init__(self, local_problem_object):
        """
            Steps of local upscaling:
        -> Assembly of local problems OK
        -> Apply boundary conditions to local problems OK
        -> Solve local problems OK
        -> Calculate upscaled permeability or upscaled transmissibility
        -> Assembly of global coarse problem
        -> Apply boundary conditions to global coarse problem
        -> Solve global coarse problem
        """
        lp = local_problem_object(preprocessor = impress, coarse_config = coarse_config, mesh_file = 'mesh/25.h5m', boundary_condition_type = 1)

    def upscaled_permeability(self):
        pass

    def upscaled_transmissibility(self):
        pass

    def assembly(self):
        pass

    def boundary_conditions(self):
        # Read data
        pass

    def solver(self):
        pass
