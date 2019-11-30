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

        print('Local uscaling class initialized')

        # Preprocessing mesh with IMPRESS
        start = time.time()
        self.mesh = preprocessor(mesh_file, dim = 3) # IMPRESS' object
        end = time.time()
        print("\nThe preprocessing step lasted {0}s".format(end-start))

        # Setting variables
        self.boundary_condition_type = boundary_condition_type
        self.coarse = self.mesh.coarse

        # Coarsening informations
        self.coarse_config = coarse_config()
        self.tree = self.coarse_config.tree['Simple']
        self.get_coarse_informations()
        self.get_coarse_centers()

        # Setting boundary conditions
        #self.set_boundary_conditions('x')

        # Coordinate System
        self.x = np.array([1,0,0]);
        self.y = np.array([0,1,0]);
        self.z = np.array([0,0,1]);

    def set_boundary_conditions(self, direction):
        """
            Indicates which function must be executed to set boundary condition on the mesh acording to the option informed.
        """
        self.direction = direction
        self.boundary_conditions = {
            1: 'self.fixed_constant_pressure()', # Fixed constant pressure
            2: 'self.fixed_linear_pressure()',   # Fixed linear pressure
            3: 'self.periodic_pressure()'        # Periodic pressure
            }

        self.bc = self.boundary_conditions.get(self.boundary_condition_type, "\nprint('Invalid boundary condition')")

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

        print('\nFixed linear pressure boundary condition applied')

    def periodic_pressure(self):
        """
        Must return a mesh with boundary conditions set in the specific given volume
        """

        print('\nPeriodic pressure boundary condition applied')

    def get_coarse_informations(self):
        """
        Access coarsening informations given to IMPRESS
        """
        self.nx = self.tree['nx']
        self.ny = self.tree['ny']
        self.nz = self.tree['nz']

        print('\nCoarsening informations accessed')

    def get_coarse_centers(self):
        """
        Calculates the center of a coarse volume and stores it in a IMPRESS variable coarse_center
        """
        for i in range(len(self.coarse.elements)):
            coords = self.coarse.elements[i].nodes.coords[:]
            x_coords = coords[:,0]
            y_coords = coords[:,1]
            z_coords = coords[:,2]

            xmax = np.amax(x_coords)
            ymax = np.amax(y_coords)
            zmax = np.amax(z_coords)

            xmin = np.amin(x_coords)
            ymin = np.amin(y_coords)
            zmin = np.amin(z_coords)

            x_center = xmin + (xmax-xmin)/2
            y_center = ymin + (ymax-ymin)/2
            z_center = zmin + (zmax-zmin)/2

            self.coarse.elements[i].coarse_center[:] = np.array([x_center, y_center, z_center])

        print('\nCoarse volume centers calculated')
