"""
Implementation of features to handle refinement of well regions. Actually, the goal of this module is to enrure that well regions will continue just as in the original fine model of reservoir. It is assumed that the vertices of the region to be refined are
coincidents with the limits of the coarse volumes.
"""
import numpy as np
from imex_integration.refinement import MeshRefinement

class UpscalingRefinement(MeshRefinement):

    def identify_fine_volumes_in_refinement_regions(self):

        volumes_in_regions = []
        self.read_refinement_info()
        centers = self.mesh.volumes.center[:]

        for region in self.refinements:
            fine_volumes = []
            fine_volumes_in_region = []
            region_coords, origin, last = self.get_refinement_coords(region)

            for direction in np.arange(3):
                greatest = np.amax(region_coords[:, direction])
                smallest = np.amin(region_coords[:, direction])
                greater = np.where(centers[:,direction] < greatest)[0]
                smaller = np.where(centers[:,direction] > smallest)[0]
                intersection = np.intersect1d(greater, smaller)
                fine_volumes.append(intersection)

            fine_volumes_in_region = np.intersect1d(fine_volumes[0], fine_volumes[1])
            fine_volumes_in_region = np.intersect1d(fine_volumes_in_region, fine_volumes[2])
            volumes_in_regions.append(fine_volumes_in_region)

        return volumes_in_regions

    def identify_coarse_volumes_in_refinement_regions(self):

        dont_upscale = np.array([], dtype = int)
        volumes_in_regions = self.identify_fine_volumes_in_refinement_regions()

        for i in range(len(self.coarse.elements)):
            ids = self.coarse.elements[i].volumes.father_id[:]
            for region in volumes_in_regions:
                intersection = np.intersect1d(ids, region)
                if intersection.size is 0:
                    pass
                else:
                    dont_upscale = np.append(dont_upscale, i)
                    break

        return dont_upscale
