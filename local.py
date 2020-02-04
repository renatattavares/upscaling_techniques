#from impress.preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as impress
from upscaling_procedures.local_upscaling import LocalUpscaling
#from imex_integration.read_dataset import read_dataset

#perm = read_dataset('imex_datasets/super.dat')
# = impress('generated_mesh.h5m', dim = 3)
lu = LocalUpscaling(mesh_file = 'mesh/25.h5m', boundary_condition_type = 1)
