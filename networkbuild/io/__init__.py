# -*- coding: utf-8 -*-

import ogr
import numpy as np
import networkx as nx
import pandas as pd
from networkbuild.classes import GeoGraph

"""
Package for reading/writing networkx based GeoGraphs
Note:  these wrap existing networkx functions for custom behavior
"""

def load_shp(shp_path):
    """ loads a shapefile into a networkx based GeoGraph object

    Parameters
    ----------

    shp_path:  string path to a line or point shapefile

    Returns
    -------

    geograph:  GeoGraph 

    """

    driver = ogr.GetDriverByName('ESRI Shapefile')
    shp = driver.Open(shp_path)
    layer = shp.GetLayer()
    spatial_ref = layer.GetSpatialRef()
    proj4 = spatialRef.ExportToProj4()

    g = nx.read_shp(shp_path)
    coords = np.array(g.nodes())

    g = nx.convert_node_labels_to_integers(g)

    return GeoGraph(srs=proj4, coords=coords, data=g)


def write_shp(geograph, shp_dir):
    """ writes a shapefile from a networkx based GeoGraph object

    Parameters
    ----------

    geograph: GeoGraph object
    shp_dir:  string path to dir to write shape files

    """

    assert geograph.node.keys() == range(geograph.coords),\
           "GeoGraph nodes and coords are not aligned"

    # need to relabel nodes by their coords
    coord_tups = map(tuple, geograph.coords)
    tup_map = {i: tup for i, tup in enumerate(coord_tups)}
    
    nx_coord_graph = nx.relabel_nodes(geograph, tup_map)

    nx.write(nx_coord_graph, shp_dir)



