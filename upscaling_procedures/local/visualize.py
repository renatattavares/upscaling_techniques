import yaml
import numpy as np
from imex_integration.mesh_constructor import MeshConstructor
from impress.preprocessor.meshHandle.finescaleMesh import FineScaleMesh

class Visualize:

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

        por, perm = self.export_info()

        coarse_model.kefx[:] =  perm[:,0]
        coarse_model.kefy[:] =  perm[:,1]
        coarse_model.kefz[:] =  perm[:,2]
        coarse_model.poref[:] =  perm[:,0]

        coarse_model.core.print()

    #
    # def print_fine_model(self):
    #
    #     perm = self.mesh.permeability[:]
    #
    #     self.mesh.permx[:] = perm[:,0]
    #     self.mesh.permy[:] = perm[:,1]
    #     self.mesh.permz[:] = perm[:,2]
    #     self.mesh.core.print()

    def export_info(self):

        effective_permeability = np.zeros((len(self.mesh.coarse.elements), 3))
        effective_porosity = np.zeros((len(self.mesh.coarse.elements), 1))

        for infos, volumes in zip(self.info, self.distribution):
            for info, volume in zip(infos, volumes):
                    effective_permeability[volume][0] = info[0]
                    effective_permeability[volume][1] = info[1]
                    effective_permeability[volume][2] = info[2]
                    effective_porosity[volume] = info[3]

        return effective_porosity, effective_permeability

    def save_info(self):

        self.effective_porosity, self.effective_permeability = self.export_info()
