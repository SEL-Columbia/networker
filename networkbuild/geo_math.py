# -*- coding: utf-8 -*-

import numpy as np
from numba import jit

"""
Module for geometric/geographic utility functions
"""

MEAN_EARTH_RADIUS_M = 6371010

def spherical_projection(coords, radius=MEAN_EARTH_RADIUS_M):
    """
    Transform angular coordinates on a sphere to x, y, z coordinates
    with origin at the center of the sphere

    Parameters:
    coords:  nx2 array of angular (degrees) longitude, latitudes
    radius:  the radius of the sphere on which to project

    Returns:
    array:   nx3 array of cartesian x, y, z coordinates with origin
        at the center of the earth (depending on the radius)

      (z) |  /
          | /.(λ,φ)
          |/ |  
     -----+-------
         /|\ |  (y) 
        / | \|
    (x)/  |  (h) 

    """
   
    assert np.shape(coords)[-1] == 2, "coords last dim must be 2 (lon, lat)"
    # transpose to operate on easily and xform to radians
    lon_rad, lat_rad = np.transpose(coords) * (np.pi / 180.0)

    # the adjacent side of the Δ between x,y plane and x,y,z point
    # (this is the hypotenuse of the Δ on x,y plane defined by λ)
    h = np.cos(lat_rad)

    # use h to calc x, y points based on λ
    x = h * np.cos(lon_rad) 
    y = h * np.sin(lon_rad)
    z = np.sin(lat_rad)

    # transpose to nx3 and project from unit circle to sphere via radius
    return np.array([x, y, z]).T * radius


def spherical_distance(coord_pairs, radius=MEAN_EARTH_RADIUS_M):
    """
    Calculate distance on sphere between pairs of coordinates

    Parameters:
    coord_pairs:  nx2x2 array of angular (degrees) longitude, latitude pairs
       (i.e. each row has a pair of coordinates to calculate dist between)
    radius:  the radius of the sphere on which to project

    Returns:
    array:   array of distances in same units as radius

    """

    assert np.shape(coord_pairs)[1:] == (2, 2), "coord_pairs must be nx2x2"
    # convert to 3-space 
    # make radius=1 because we want unit vectors to make below math easier
    coords_3 = spherical_projection(coord_pairs, radius=1)

    # law of cosines and vector algebra:
    # cosθ = u·v/(||u||*||v||)
    # in this case, ||u|| and ||v|| are 1 (on the unit sphere)
    cos_theta = np.sum(np.product(coords_3, axis=1), axis=1)
    theta = np.arccos(cos_theta)
    return radius * theta

def spherical_distance_scalar(coord_pair, radius=MEAN_EARTH_RADIUS_M):
    """
    Wrapper for spherical_distance which takes a single set of pairs

    Parameters:
    coord_pair:  1x2 array of lon, lat pairs (e.g. [[30, 30], [31, 31]])
    radius:  radius of sphere on which to project

    Returns:
    distance:  single distance between the pair
    """
    
    return spherical_distance([coord_pair], radius)[0]

