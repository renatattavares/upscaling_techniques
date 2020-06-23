from imex_integration.read_dataset import read_dataset
from upscaling_procedures.local.local_problems import LocalProblems
from upscaling_procedures.local.local_upscaling import LocalUpscaling
from upscaling_procedures.local.parallel_local_problems import ParallelLocalProblems
from upscaling_procedures.local.parallel_local_upscaling import ParallelLocalUpscaling
from impress.preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as impress

############### RUN LOCAL PROBLEMS ###############
#lp = LocalProblems(mesh_file = 'mesh/20.h5m', dataset = None)
#lp = ParallelLocalProblems(mesh_file = 'mesh/20.h5m', dataset = None)
#lp = LocalProblems(mesh_file = None, dataset = 'imex_datasets/data_20_mesh.dat')

############### RUN LOCAL UPSCALING ###############
#lu = LocalUpscaling(mesh_file = None, dataset = 'imex_datasets/data_20_mesh.dat')
lu = LocalUpscaling(mesh_file = 'mesh/20.h5m', dataset = None)

############### RUN PARALLEL LOCAL UPSCALING ###############
#lu = ParallelLocalUpscaling(mesh_file = 'mesh/20.h5m', dataset = None)
#lu = ParallelLocalUpscaling(mesh_file = None, dataset ='imex_datasets/spe_10_case_2/spe_10_case_2.dat')

############### READ DATASET ###############
#porosity, permeability = read_dataset('imex_datasets/spe_10_case_2/spe_10_case_2.dat')

############### WRITE DATASET ###############
#new_dataset = DatasetWriter('imex_datasets/spe_10_case_2/spe_10_case_2.dat')

############### RUN IMPRESS ###############
#mesh = impress(mesh_file = 'generated_mesh.h5m', dim = 3)
