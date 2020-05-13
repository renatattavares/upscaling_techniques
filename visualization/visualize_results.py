import yaml
from imex_integration.mesh_constructor import MeshConstructor as MC

def create_coarse_mesh(file_path):

    with open(file_path, 'r') as file:
        data = yaml.safe_load(file)

    dict = data['Simple']
    nx = dict['nx']
    ny = dict['ny']
    nz = dict['nz']

    number_elements = [nx, ny, nz]
    length_elements = [1, 1 ,1]

    generator = MC(number_elements, length_elements)

def print_results():
    pass
