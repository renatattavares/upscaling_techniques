import unittest
import numpy as np
from upscaling_procedures.local.parallel_local_upscaling import ParallelLocalUpscaling

class Test20Case(unittest.TestCase):

    def test_effective_permeability(self):

        self.lu = ParallelLocalUpscaling(mesh_file = 'mesh/20.h5m', dataset = None)
        self.lu.save_info()

        self.assertEqual(self.lu.effective_permeability[:,0].all(), np.ones(125).all())
        self.assertEqual(self.lu.effective_permeability[:,1].all(), np.ones(125).all())
        self.assertEqual(self.lu.effective_permeability[:,2].all(), np.ones(125).all())
