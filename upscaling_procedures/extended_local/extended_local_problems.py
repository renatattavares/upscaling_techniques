"""
Module for treatment of extended local problems to apply local upscaling technique in structured tridimensional meshes
"""
from upscaling_procedures.local.local_problems import LocalProblems

class ExtendedLocalProblems(LocalProblems):

    def __init__(self, mesh_file = None, dataset = None):
        initial_time = time.time()

        print('\n##### ExtendedLocalProblems class initialized #####')

        super().__init__(mesh_file, dataset)

        
