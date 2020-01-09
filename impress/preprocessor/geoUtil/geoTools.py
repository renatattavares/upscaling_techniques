"""
Geometric methods to compute volumes, areas, distances of the mesh entities
"""
# Geoemtric Module
# Create by Artur Castiel and Renata Tavares

#import pdb
import numpy as np

def normal_vec_2d(coords0,coords1):
    vec = coords1 - coords0
    norm = np.linalg.norm(vec, axis = 1)
    norm = 1/norm
    return np.array([vec[:,1], -vec[:,0], vec[:,2] ]).T * norm[:,np.newaxis]
    # distance = (np.inner(vec, vec, axis = 0))

def normal_vec(coords0, coords1, coords2):
    vec1 = coords1 - coords0
    vec2 = coords2 - coords0
    #cross = cross_numba(vec1,vec2)
    cross = np.cross(vec1,vec2)
    norm = np.linalg.norm(cross, axis = 1)
    norm = 1/norm
    return  cross * norm[:,np.newaxis]
    # a = cross * norm

def point_distance(coords_1, coords_2):
    dist_vector = coords_1 - coords_2
    distance = sqrt(np.dot(dist_vector, dist_vector))
    return distance

def get_average(coords_list):
    N = len(coords_list)
    return sum(coords_list)*(1/N)
