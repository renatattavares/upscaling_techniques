from imex_integration.read_dataset import read_dataset
from upscaling_procedures.local.local_problems import LocalProblems
from upscaling_procedures.local.local_upscaling import LocalUpscaling
from upscaling_procedures.local.parallel_local_problems import ParallelLocalProblems
from upscaling_procedures.local.parallel_local_upscaling import ParallelLocalUpscaling
from impress.preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as impress

lu = ParallelLocalUpscaling(mesh_file = None, dataset ='imex_datasets/spe_10_case_2/spe_10_case_2.dat')
