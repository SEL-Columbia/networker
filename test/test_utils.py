# -*- coding: utf-8 -*-

import numpy as np
import networker.utils as utils

def test_array2d_coords_dict_conversions():
    """
    Test whether we convert 2d arrays to dicts and back OK
    """

    coords_dict = {'g1': (1.0, 1.0), 
                   'g2': (2.0, 2.0),
                   'g3': (3.0, 3.0)}


    test_array2d, index_map = utils.coords_dict_to_array2d(coords_dict)

    test_coords_dict = utils.array2d_to_coords_dict(test_array2d, index_map)
    assert test_coords_dict == coords_dict, \
        "array to coords conversion error"
