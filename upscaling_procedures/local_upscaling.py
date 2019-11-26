"""
Module of local upscaling technique in structured tridimensional meshes
"""

import numpy as np

class local_upscaling:

    def __init__(self, mesh_file, mesh = None):
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

    self.mesh = mesh(mesh_file, dim = 3) # IMPRESS' object

    def set_boundary_conditions(self, boundary_condition_type):
        """
        Indicates which function must be executed to set boundary condition on the mesh acording to the option informed.
        """

        self.boundary_conditions = {
            'fixed constant pressure': 'self.fixed_constant_pressure()',
            'fixed linear pressure': 'self.fixed_linear_pressure()',
            'periodic pressure': 'self.periodic_pressure()'
            }

        self.bc = self.boundary_conditions.get(boundary_condition_type, "print('Invalid boundary condition')")

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
        print('Fixed linear pressure boundary condition applied')

    def periodic_pressure(self):
        print('Periodic pressure boundary condition applied')
