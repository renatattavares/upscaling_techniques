import yaml
import numpy as np
from pymoab import core, types, rng
from imex_integration.refinement import Refinement

class MeshConstructor(Refinement):

    def __init__(self, number_elements, length_elements, mesh_file):

        print('\n##### Generating mesh file #####')
        self.mesh_file = mesh_file
        self.mbcore = core.Core() # MOAB Core -> Mesh Management
        refine = self.check_if_refinement_is_required()

        if refine is False: # Refinement is not required
            print('\nCreating vertices coordinates...')
            self.coords = self.create_vertices_coords(number_elements, length_elements, np.array([0,0,0]))
            print('Creating mesh connectivities...')
            self.mesh_connectivity = self.create_mesh_connectivity(number_elements, length_elements, self.coords) # Indexes of vertices coords that composes an element
            print("Creating elements' handles...")
            self.elements_handles = self.create_elements_handles(self.mesh_connectivity, self.coords)
            print('Writing file...')
            self.mbcore.write_file(self.mesh_file)
            print('\n##### Mesh file created #####')

        else: # Refinement is required
            print('\nCreating vertices coordinates...')
            self.coords = self.create_vertices_coords(number_elements, length_elements, np.array([0,0,0]))
            print('Creating mesh connectivities...')
            self.mesh_connectivity = self.create_mesh_connectivity(number_elements, length_elements, self.coords) # Indexes of vertices coords that composes an element
            self.read_refinement_info()
            print('Starting refinement step...\n')
            self.refine_regions()
            coords = np.reshape(self.coords, newshape = (int(len(self.coords)/3), 3))
            print('Rewriting mesh connectivity...')
            new_mesh_connectivity = self.rewrite_mesh_connectivity(self.mesh_connectivity, coords, self.new_coords)
            print("Creating elements' handles...")
            self.handles = self.create_elements_handles(new_mesh_connectivity, self.new_coords.flatten())
            print('Writing file...')
            self.mbcore.write_file(self.mesh_file)
            print('\n##### Mesh file created #####')

    def create_vertices_coords(self, number_elements, length_elements, mesh_origin):

        nx = number_elements[0]
        ny = number_elements[1]
        nz = number_elements[2]

        dx = length_elements[0]
        dy = length_elements[1]
        dz = length_elements[2]

        num_vertex = int((nx+1)*(ny+1)*(nz+1))
        vertex_coords = np.zeros(num_vertex*3)

        for i in range(num_vertex):
            vertex_coords[3*i] = (i % (nx+1))*dx + mesh_origin[0]
            vertex_coords[3*i+1] = ((i // (nx+1)) % (ny+1))*dy + mesh_origin[1]
            vertex_coords[3*i+2] = ((i // ((nx+1)*(ny+1))) % (nz+1))*dz + mesh_origin[2]

        return vertex_coords

    def create_mesh_connectivity(self, number_elements, length_elements, vertex_coords):

        k = 0
        nx = number_elements[0]
        ny = number_elements[1]
        nz = number_elements[2]
        dx = length_elements[0]
        dy = length_elements[1]
        dz = length_elements[2]
        num_elements = nx*ny*nz
        mesh_connectivity = np.zeros((num_elements, 8), dtype = int)

        # Defining the vertices which can start an element, i.e., starting from this
        # vertex one can get five more valid vertices to build an element.
        coords_indexes = np.arange(int(len(vertex_coords)/3), dtype = int)
        indexes = [int(v) for v in coords_indexes if (vertex_coords[(3*v)] != nx*dx) and (vertex_coords[(3*v)+1] != ny*dy) and (vertex_coords[(3*v)+2] != nz*dz)]

        for i in indexes:
            mesh_connectivity[k] = [coords_indexes[i], coords_indexes[i+1], coords_indexes[int(i+nx+2)], coords_indexes[int(i+nx+1)], coords_indexes[int(i+(nx+1)*(ny+1))], coords_indexes[int(i+(nx+1)*(ny+1)+1)],                 coords_indexes[int(i+(nx+1)*(ny+2)+1)], coords_indexes[int(i+(nx+1)*(ny+2))]]
            k += 1

        return mesh_connectivity

    def create_elements_handles(self, connectivity, vertex_coords):

        self.vertex_handles = self.mbcore.create_vertices(vertex_coords) # Criando os handles dos vértices
        connectivity = (connectivity + 1).astype('uint64')
        elements_handles = rng.Range([self.mbcore.create_element(types.MBHEX, x) for x in connectivity])

        return elements_handles
