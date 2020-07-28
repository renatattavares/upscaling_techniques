import yaml
import numpy as np
from pymoab import core, types, rng

class MeshRefinement:

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

    def refine_regions(self):

        coords = np.reshape(self.coords, newshape = (int(len(self.coords)/3), 3))
        region_number = np.arange(len(self.refinements) + 1)
        self.new_coords = coords
        refinements_connectivities = []
        refinements_coords = []

        for region, number in zip(self.refinements, region_number):
            print('Refining region number {}'.format(number))
            region_coords, origin, last  = self.get_refinement_coords(region)
            length = np.array(region['length'])
            elements = ((last - origin) / length).astype('int')
            refinement_coords = self.create_vertices_coords(elements, length, np.array([0,0,0]))
            refinements_connectivities.append(self.create_mesh_connectivity(elements, length, refinement_coords))
            refinement_coords = np.reshape(refinement_coords, newshape = (int(len(refinement_coords)/3), 3))
            refinement_coords = refinement_coords + origin
            refinements_coords.append(refinement_coords)
            self.update_mesh_connectivity(coords, refinement_coords)

        print('\nUpdating vertices coordinates...')
        for c in refinements_coords:
            coords = self.update_mesh_coords(coords, c)
            coords = np.append(coords, c.flatten())
            coords = np.reshape(coords, newshape = (int(len(coords)/3), 3))

        self.new_coords = coords

        print("Creating refinement elements' handles...")
        for connectivity, coords in zip(refinements_connectivities, refinements_coords):
            new_mesh_connectivity = self.rewrite_mesh_connectivity(connectivity, coords, self.new_coords)
            self.create_elements_handles(new_mesh_connectivity, self.new_coords.flatten())

    def get_refinement_coords(self, dict):

        region_coords = []

        for i in range(8):
            region_coords.append(dict[i+1])

        region_coords = np.array(region_coords)
        norms = np.linalg.norm(region_coords, axis = 1)
        origin = region_coords[np.argmin(norms)] # The closest point to mesh origin in refinement coords
        last = region_coords[np.argmax(norms)] # The farthest point to mesh origin in refinement coords

        return region_coords, origin, last

    def update_mesh_coords(self, coords, refinement_coords):

        delete = []
        new_coords = coords

        for i in refinement_coords:
            k = 0
            for j in coords:
                if np.array_equal(i,j):
                    delete.append(k)
                k += 1

        new_coords = np.delete(coords, delete, axis = 0)

        return new_coords

    def update_mesh_connectivity(self, coords, refinement_coords):

        delete = []
        for element in self.mesh_connectivity:
            #print('Elemento')
            vertices = coords[element]
            d = True
            for v in vertices:
                array = np.full_like(refinement_coords, v)
                d = np.any(np.all(array == refinement_coords, axis = 1))
                if d == False:
                    break
            if d == True:
                delete.append(np.where(self.mesh_connectivity == element)[0])

        self.mesh_connectivity = np.delete(self.mesh_connectivity, delete, axis = 0)

    def rewrite_mesh_connectivity(self, mesh_connectivity, old_coords, new_coords):

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
