"""
Module for treatment of local problems to apply local upscaling technique in structured tridimensional meshes
"""
import time
import numpy as np
#import xlsxwriter
from local_problems.solver_local_problems import Solver
from local_problems.boundary_conditions import BoundaryConditions
from local_problems.assembly import Assembly
from impress.preprocessor.meshHandle.configTools.configClass import coarseningInit as coarse_config
from impress.preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as impress

class LocalProblems(BoundaryConditions, Solver, Assembly):
    @profile
    def __init__(self, mesh_file = None, boundary_condition_type = None):

        print('\n##### Treatment of local problems #####')

        # Preprocessing mesh with IMPRESS
        print('\nPre-processing mesh with IMPRESS...')
        start = time.time()
        self.mesh = impress(mesh_file, dim = 3) # IMPRESS' object
        end = time.time()
        print("\nThis step lasted {0}s".format(end-start))

        # Setting variables
        print('\nAccessing coarsening informations from IMPRESS and setting important variables...')

        self.mesh.permeability[:] = np.array([1,1,1])
        self.get_mesh_informations(coarse_config())
        self.coarse = self.mesh.coarse
        self.boundary_condition_type = boundary_condition_type
        self.number_coarse_volumes = len(self.coarse.elements) # Number of volumes from the coarse mesh
        self.number_volumes_local_problem = len(self.mesh.volumes)/(self.nx*self.ny*self.nz) # Number of fine scale volumes inside a coarse volume
        self.x = np.array([1,0,0])
        self.y = np.array([0,1,0])
        self.z = np.array([0,0,1])
        self.pressure_gradient = 1
        self.direction_string = np.array(['x', 'y', 'z'])
        self.directions_dictionary = {
            'x': self.x,
            'y': self.y,
            'z': self.z
            }
        # Set and solve local problems in x, y and z directions
        for i in self.direction_string:
            print('\nAssembly of local problems in {} direction...'.format(i))
            start = time.time()
            self.assembly_local_problem()
            end = time.time()
            print("\nThis step lasted {}".format(end-start))

            print('\nSetting boundary conditions in {} direction...'.format(i))
            start = time.time()
            self.set_boundary_conditions(i)
            end = time.time()
            print("\nThis step lasted {}".format(end-start))

            print('\nSolving local problems in {} direction...'.format(i))
            start = time.time()
            self.solve_local_problems()
            end = time.time()
            print("\nThis step lasted {}".format(end-start))

            #return self.pressure_x, self.pressure_y, self.pressure_z
