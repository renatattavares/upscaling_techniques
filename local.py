import numpy as np
from imex_integration.read_dataset import read_dataset
from imex_integration.write_dataset import DatasetWriter
from upscaling_procedures.local.local_problems import LocalProblems
from upscaling_procedures.local.local_upscaling import LocalUpscaling
from upscaling_procedures.local.parallel_local_upscaling import ParallelLocalUpscaling
from upscaling_procedures.extended_local.extended_local_problems import ExtendedLocalProblems
from upscaling_procedures.extended_local.extended_local_upscaling import ExtendedLocalUpscaling
from upscaling_procedures.extended_local.parallel_extended_local_upscaling import ParallelExtendedLocalUpscaling
from impress.preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as impress

#from impress.preprocessor.meshHandle.finescaleMesh import FineScaleMesh as impress

############### RUN LOCAL PROBLEMS ###############
#lp = LocalProblems(mesh_file = 'mesh/20.h5m', dataset = None)
#lp = LocalProblems(mesh_file = None, dataset = 'imex_datasets/data_20_mesh.dat')

############### RUN LOCAL UPSCALING ###############
#lu = LocalUpscaling(mesh_file = None, dataset = 'imex_datasets/data_20_mesh.dat')
#lu = LocalUpscaling(mesh_file = 'mesh/20.h5m', dataset = None)

############### RUN PARALLEL LOCAL UPSCALING ###############
#lu = ParallelLocalUpscaling(mesh_file = 'mesh/20.h5m', dataset = None)
#lu = ParallelLocalUpscaling(mesh_file = None, dataset ='imex_datasets/spe_10_case_2/spe_10_case_2.dat')
#por, perm = lu.export_info()
#info_refined_volumes = lu.export_info_refined_volumes()

############### RUN PARALLEL EXTENDED LOCAL UPSCALING ###############
#elu = ExtendedLocalProblems(mesh_file = 'mesh/20.h5m', dataset = None)
#elu = ExtendedLocalUpscaling(mesh_file = 'mesh/20.h5m', dataset = None)
elu = ParallelExtendedLocalUpscaling(mesh_file = 'mesh/20.h5m', dataset = None)
#elu = ParallelExtendedLocalUpscaling(mesh_file = None, dataset = 'imex_datasets/spe_10_case_2/spe_10_case_2.dat')
#por, perm = elu.export_info()
#info_refined_volumes = elu.export_info_refined_volumes()

############## READ DATASET ###############
# por, perm, number_elements, length_elements  = read_dataset('imex_datasets/spe_10_case_2/spe_10_case_2.dat', 'dataset_mesh.h5m')
# info_refined_volumes = np.array([1,2,3,4])
#
# ############### WRITE DATASET ###############
# new_dataset = DatasetWriter('imex_datasets/spe_10_case_2/spe_10_case_2.dat', np.array([12,44,17]), np.array([60,220,85]), np.array([20,10,2]), por, perm, info_refined_volumes)
#new_dataset = DatasetWriter('imex_datasets/spe_10_case_2/spe_10_case_2.dat', lu.coarsening, lu.number_elements, lu.length_elements, por, perm, info_refined_volumes)

############### RUN IMPRESS ###############
#mesh = impress(mesh_file = 'mesh/20.h5m', dim = 3)
