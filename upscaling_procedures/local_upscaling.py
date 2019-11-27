"""
Module of local upscaling technique in structured tridimensional meshes
"""
import time
import numpy as np

class LocalUpscaling:

    def __init__(self, preprocessor, coarse_config, mesh_file = None, boundary_condition_type = None):
        """
            Steps of local upscaling:

        -> Apply boundary conditions to local problems
        -> Assembly of local problems
        -> Solve local problems
        -> Calculate upscaled permeability or upscaled transmissibility
        -> Apply boundary conditions to global coarse problem
        -> Assembly of global coarse problem
        -> Solve global coarse problem
        """

        # Preprocessing mesh with IMPRESS
        start = time.time()
        self.mesh = preprocessor(mesh_file, dim = 3) # IMPRESS' object
        end = time.time()
        print("The preprocessing step lasted {0}s".format(end-start))

        # Coarsening informations
        self.coarse_config = coarse_config()
        self.tree = self.coarse_config.tree['Simple']
        self.get_coarse_informations()
        print('Coarsening informations accessed')

        # Setting class variables
        self.boundary_condition_type = boundary_condition_type


    def set_boundary_conditions(self):
        """
            Indicates which function must be executed to set boundary condition on the mesh acording to the option informed.
        """
        self.boundary_conditions = {
            1: 'self.fixed_constant_pressure()', # Fixed constant pressure
            2: 'self.fixed_linear_pressure()',   # Fixed linear pressure
            3: 'self.periodic_pressure()'        # Periodic pressure
            }

        self.bc = self.boundary_conditions.get(self.boundary_condition_type, "print('Invalid boundary condition')")

        exec(self.bc)

    def assembly(self):
        pass

    def solver(self):
        pass

    def upscaled_permeability(self):
        pass

    def upscaled_transmissibility(self):
        pass

    def fixed_constant_pressure(self):
        """
        Must return a mesh with boundary conditions set in the specific given volume
        """
        
        print('Fixed constant pressure boundary condition applied')

    def fixed_linear_pressure(self):
        """
        Must return a mesh with boundary conditions set in the specific given volume
        """
        print('Fixed linear pressure boundary condition applied')

    def periodic_pressure(self):
        """
        Must return a mesh with boundary conditions set in the specific given volume
        """
        print('Periodic pressure boundary condition applied')

    def get_coarse_informations(self):
        """
        Access coarsening informations given to IMPRESS
        """
        self.nx = self.tree['nx']
        self.ny = self.tree['ny']
        self.nz = self.tree['nz']
