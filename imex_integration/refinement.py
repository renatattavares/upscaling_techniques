import yaml
import numpy as np
from pymoab import core, types, rng

class Refinement():

    def check_if_refinement_is_required(self, refinement_file = 'refinement.yml'):

        with open(refinement_file, 'r') as file:
            self.data = yaml.safe_load(file)

        refine = self.data['refine']

        return refine

    def read_refinement_info(self):

        self.refinements = []

        for key in self.data.keys():
            if key == 'refine':
                pass
            else:
                self.refinements.append(self.data[key])

    def create_refinement_vertices_coords(self):

        coords = np.reshape(self.coords, newshape = (int(len(self.coords)/3), 3))
        region_coords = []
        self.connectivities = []

        for region in self.refinements:
            delete = []

            for i in range(8):
                region_coords.append(region[i+1])

            length = np.array(region['length'])
            elements = ((np.array(region_coords[7]) - np.array(region_coords[0])) / length).astype('int')
            refinement_coords = self.create_vertices_coords(elements, length, region_coords[0])
            self.refinement_coords = refinement_coords
            refinement_connectivity = self.create_mesh_connectivity(elements, length, refinement_coords)
            refinement_coords = np.reshape(refinement_coords, newshape = (int(len(refinement_coords)/3), 3))

            for element in self.mesh_connectivity:
                vertices = coords[element]
                delete_element = np.all(np.isin(vertices, refinement_coords))
                if delete_element == True:
                    delete.append(np.where(self.mesh_connectivity == element)[0])

            self.mesh_connectivity = np.delete(self.mesh_connectivity, delete, axis = 0)
            self.connectivities.append(refinement_connectivity)

        region_coords = np.array(region_coords)
        delete = []
        for i in region_coords:
            k = 0
            for j in coords:
                if np.array_equal(i,j):
                    delete.append(k)
                k += 1

        self.new_coords = np.delete(coords, delete, axis = 0).flatten()
        self.new_coords = np.append(self.new_coords, refinement_coords.flatten())


    def rewrite_mesh_connectivity(self, mesh_connectivity, old_coords, new_coords):

        old_coords = np.reshape(old_coords, newshape = (int(len(old_coords)/3), 3))
        new_coords = np.reshape(new_coords, newshape = (int(len(new_coords)/3), 3))
        new_mesh_connectivity = []

        for element in mesh_connectivity:
            element_connectivity = []
            vertices = old_coords[element]
            for vertice in vertices:
                i = 0
                for coord in new_coords:
                    if np.array_equal(vertice, coord):
                        element_connectivity.append(i)
                        break
                    i += 1
            new_mesh_connectivity.append(element_connectivity)

        return np.array(new_mesh_connectivity)
