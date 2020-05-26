import yaml
import numpy as np
from imex_integration.mesh_constructor import MeshConstructor
from impress.preprocessor.meshHandle.finescaleMesh import FineScaleMesh

class Visualize:

    def __init__(self):
        pass

    def print_coarse_model(self, mesh_file = 'mesh/coarse_model.h5m'):

        fine_number_elements = self.number_elements
        fine_length_elements = self.length_elements
        self.coarse_number_elements = np.array([], dtype = int)

        with open('impress/input_cards/coarsening.yml', 'r') as coarsening:
            new_mesh_data = yaml.safe_load(coarsening)

        coarsening = new_mesh_data['Simple']

        for key in coarsening.keys():
            self.coarse_number_elements = np.append(self.coarse_number_elements,coarsening[key])

        self.coarse_length_elements = (fine_number_elements/self.coarse_number_elements)*fine_length_elements
        generator = MeshConstructor(self.coarse_number_elements, self.coarse_length_elements, mesh_file)
        coarse_model = FineScaleMesh(mesh_file)

        for a, b in zip(self.effective_permeability, self.distribution):
            for c, d in zip(a, b):
                coarse_model.kefx[int(d)] = c[0]
                coarse_model.kefy[int(d)] = c[1]
                coarse_model.kefz[int(d)] = c[2]

        coarse_model.core.print()

    def print_fine_model(self):

        perm = self.mesh.permeability[:]

        self.mesh.permx[:] = perm[:,0]
        self.mesh.permy[:] = perm[:,1]
        self.mesh.permz[:] = perm[:,2]
        self.mesh.core.print()
