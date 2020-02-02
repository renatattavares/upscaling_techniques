#from upscaling_procedures.local_upscaling import LocalUpscaling
from imex_integration.read_dataset import read_dataset

read_dataset('imex_datasets/spe10.dat')
#lu = LocalUpscaling(mesh_file = 'generated_mesh.h5m', boundary_condition_type = 1)
