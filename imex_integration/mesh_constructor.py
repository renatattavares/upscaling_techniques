import numpy as np
from pymoab import core, types, rng

class MeshConstructor:
    def __init__(self, number_elements, length_elements):
        print('\n##### Generating mesh file #####')
        elements = number_elements
        length = length_elements

        self.dx = length[0]
        self.dy = length[1]
        self.dz = length[2]

        self.nx = elements[0]
        self.ny = elements[1]
        self.nz = elements[2]

        self.mbcore = core.Core() # Criando instância da classe Core que gerencias as operações na malha
        self.num_elements = self.nx*self.ny*self.nz
        num_vertex = (self.nx+1)*(self.ny+1)*(self.nz+1)

        print('\nCreating vertices coordinates...')
        self.vertex_coords = np.zeros(num_vertex*3)
        for i in range(num_vertex):
            self.vertex_coords[3*i] = (i % (self.nx+1))*self.dx
            self.vertex_coords[3*i+1] = ((i // (self.nx+1)) % (self.ny+1))*self.dy
            self.vertex_coords[3*i+2] = ((i // ((self.nx+1)*(self.ny+1))) % (self.nz+1))*self.dz

        self.vertex_handles = self.mbcore.create_vertices(self.vertex_coords) # Criando os handles dos vértices
        print('Creating connectivities...')
        self.mesh_connectivity = self.create_mesh_connectivity()

        print('Creating element handles...')
        self.elem_handles = rng.Range([self.mbcore.create_element(types.MBHEX, x) for x in self.mesh_connectivity])
        self.write_files()

    def create_mesh_connectivity(self):

        k = 0
        self.mesh_connectivity = np.zeros((self.num_elements, 8), dtype=np.uint64)

        # Defining the vertices which can start an element, i.e., starting from this
        # vertex one can get five more valid vertices to build an element.
        indexes = [v-1 for v in self.vertex_handles
                        if (self.vertex_coords[3*(int(v)-1)] != self.nx*self.dx) and
                            (self.vertex_coords[3*(int(v)-1)+1] != self.ny*self.dy) and
                            (self.vertex_coords[3*(int(v)-1)+2] != self.nz*self.dz)]

        for i in indexes:
            self.mesh_connectivity[k] = [self.vertex_handles[int(i)], self.vertex_handles[int(i)+1], self.vertex_handles[int(i+self.nx+2)], self.vertex_handles[int(i+self.nx+1)], self.vertex_handles[int(i+(self.nx+1)*(self.ny+1))], self.vertex_handles[int(i+(self.nx+1)*(self.ny+1)+1)],                 self.vertex_handles[int(i+(self.nx+1)*(self.ny+2)+1)], self.vertex_handles[int(i+(self.nx+1)*(self.ny+2))]]

            k += 1

        return self.mesh_connectivity

    def write_files(self):
        self.meshset = self.mbcore.create_meshset()
        self.mbcore.add_entities(self.meshset, self.elem_handles)

        # Escrevendo malha em arquivo vtk para visualização no visit. Para utilizar a função write_file é necessário ter uma entidade iterável. Portanto, é necessária a criação de um range
        mesh_file = rng.Range(self.meshset)
        self.mbcore.write_file('mesh/dataset_mesh.h5m')
