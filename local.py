import numpy as np
from imex_integration.read_dataset import read_dataset
from imex_integration.write_dataset import DatasetWriter
from upscaling_procedures.local.local_problems import LocalProblems
from upscaling_procedures.local.local_upscaling import LocalUpscaling
from upscaling_procedures.local.parallel_local_upscaling import ParallelLocalUpscaling
from impress.preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as impress

#from impress.preprocessor.meshHandle.finescaleMesh import FineScaleMesh as impress

############### RUN LOCAL PROBLEMS ###############
#lp = LocalProblems(mesh_file = 'mesh/20.h5m', dataset = None)
#lp = LocalProblems(mesh_file = None, dataset = 'imex_datasets/data_20_mesh.dat')

############### RUN LOCAL UPSCALING ###############
#lu = LocalUpscaling(mesh_file = None, dataset = 'imex_datasets/data_20_mesh.dat')
#lu = LocalUpscaling(mesh_file = 'mesh/20.h5m', dataset = None)

############### RUN PARALLEL LOCAL UPSCALING ###############
lu = ParallelLocalUpscaling(mesh_file = 'mesh/20.h5m', dataset = None)
#lu = ParallelLocalUpscaling(mesh_file = None, dataset ='imex_datasets/spe_10_case_2/spe_10_case_2.dat')
#por, perm = lu.export_info()

############## READ DATASET ###############
#por, perm, number_elements, length_elements  = read_dataset('imex_datasets/spe_10_case_2/spe_10_case_2.dat', 'mesh1.h5m')


############### WRITE DATASET ###############
#new_dataset = DatasetWriter('imex_datasets/spe_10_case_2/spe_10_case_2.dat', np.array([1,1,1]), np.array([1,1,1]), por, perm)
#lenght = (lu.number_elements/lu.coarsening)*lu.length_elements
#new_dataset = DatasetWriter('imex_datasets/spe_10_case_2/spe_10_case_2.dat', lu.coarsening, lenght, por, perm)

############### RUN IMPRESS ###############
#mesh = impress(mesh_file = 'mesh/10.h5m', dim = 3)
