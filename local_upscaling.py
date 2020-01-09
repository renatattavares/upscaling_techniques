from impress.preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as impress
from impress.preprocessor.meshHandle.configTools.configClass import coarseningInit as coarse_config
from upscaling_procedures.local_upscaling import LocalUpscaling, LocalProblems

lu = LocalUpscaling(LocalProblems)
