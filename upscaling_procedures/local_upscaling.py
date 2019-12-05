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

        print('Local uscaling class initialized\n')

        # Preprocessing mesh with IMPRESS
        start = time.time()
        self.mesh = preprocessor(mesh_file, dim = 3) # IMPRESS' object
        end = time.time()
        print("\nThe preprocessing step lasted {0}s".format(end-start))

        # Setting variables
        self.boundary_condition_type = boundary_condition_type
        self.coarse = self.mesh.coarse
        self.number_coarse_volumes = len(self.coarse.elements)

        # Coordinate System
        self.x = np.array([1,0,0]);
        self.y = np.array([0,1,0]);
        self.z = np.array([0,0,1]);

        # Coarsening informations
        self.coarse_config = coarse_config()
        self.tree = self.coarse_config.tree['Simple']
        self.get_coarse_informations()
        self.get_coarse_centers()

        # Setting boundary conditions
        self.set_boundary_conditions('x', 0)

    def set_boundary_conditions(self, direction, local_problem):
        """
            Indicates which function must be executed to set boundary condition on the mesh acording to the option informed.
        """
        self.local_problem = local_problem
        self.direction = direction

        self.directions = {
            'x': self.x,
            'y': self.y,
            'z': self.z
            }

        self.boundary_conditions = {
            1: 'self.fixed_constant_pressure()', # Fixed constant pressure
            2: 'self.fixed_linear_pressure()',   # Fixed linear pressure
            3: 'self.periodic_pressure()'        # Periodic pressure
            }

        self.direction = self.directions.get(self.direction)

        self.bc = self.boundary_conditions.get(self.boundary_condition_type, "\nprint('Invalid boundary condition')")

        self.get_coarse_neighbors()

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
        self.neighbors = np.amax(self.correct_neighbors, axis = 1)

        coarse_neighbors = self.coarse.iface_neighbors(self.local_problem)[0]
        coarse_faces = self.coarse.iface_neighbors(self.local_problem)[1]
        boundary_element = np.isin(coarse_neighbors, self.neighbors[self.local_problem])
        index = np.where(boundary_element == True)
        interface = coarse_faces[index]
        faces = self.coarse.interfaces_faces[interface]

        adj_volumes = self.mesh.faces.bridge_adjacencies(faces,2,3)
        adj_volumes = np.concatenate(adj_volumes, axis = 0)
        adj_volumes = np.unique(adj_volumes)
        fine_volumes = self.coarse.elements[self.local_problem].volumes.father_id[:]
        correct_volumes = np.isin(fine_volumes, adj_volumes)
        index =  np.where(correct_volumes == True)
        volumes = adj_volumes[index]

        self.mesh.pressure[volumes] = 1000
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
        for i in range(self.number_coarse_volumes):
            coords = self.coarse.elements[i].nodes.coords[:]
            x_coords = coords[:,0]
            y_coords = coords[:,1]
            z_coords = coords[:,2]

            xmax = x_coords.max()
            ymax = y_coords.max()
            zmax = z_coords.max()

            xmin = x_coords.min()
            ymin = y_coords.min()
            zmin = z_coords.min()

            x_center = xmin + (xmax-xmin)/2
            y_center = ymin + (ymax-ymin)/2
            z_center = zmin + (zmax-zmin)/2

            self.coarse.elements[i].coarse_center[:] = np.array([x_center, y_center, z_center])

        print('\nCoarse volume centers calculated')

    def get_coarse_neighbors(self):

        self.correct_neighbors = np.zeros((self.number_coarse_volumes, 2))

        for i in range(self.number_coarse_volumes):
            neighbors = self.coarse.elements[i].faces.coarse_neighbors
            boundary_element = np.isin(neighbors, self.number_coarse_volumes)
            index = np.where(boundary_element == True)
            neighbors = np.delete(neighbors, index)

            for j in neighbors:
                center = self.coarse.elements[j].coarse_center[0]
                vector = self.coarse.elements[i].coarse_center[0] - center
                pv = np.cross(vector, self.direction)

                if np.linalg.norm(pv) == 0:
                    if self.correct_neighbors[i,0] == 0:
                        self.correct_neighbors[i,0] = j
                    else:
                        self.correct_neighbors[i,1] = j

        print('\nCoarse neighbors defined')
