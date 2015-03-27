import numpy as np
import networkx as nx
from nose.tools import eq_

from networkbuild.geo_math import nn_dists, ang_to_vec_coords
from networkbuild import KDTree

def simple_coord_set():

    coords = np.array([
       [ 0.35893606,  0.67506964], 
       [ 0.33262164,  0.66603501],
       [ 0.34752438,  0.68935952]])

    return coords


def TestKDTree():

    simple_coords = simple_coord_set()
    proj_coords = ang_to_vec_coords(simple_coords)
    coord_nn_dists = nn_dists(simple_coords)
    proj_nn_dists = nn_dists(proj_coords, spherical=False)

    coord_idx = np.arange(proj_coords.shape[0])
    kdt = KDTree(proj_coords)
    kd_nn_idx = [kdt.query_subset(proj_coords[i], coord_idx[coord_idx != i])[0] for i in range(len(proj_coords))]

    assert all(coord_nn_dists[1] == proj_nn_dists[1]), "spherical and cartesian adjacency do not match"
    assert all(kd_nn_idx == proj_nn_dists[1]), "kdtree and cartesian adjacency do not match"
     
