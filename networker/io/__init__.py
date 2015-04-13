# -*- coding: utf-8 -*-

import ogr, osr
import numpy as np
import networkx as nx
import pandas as pd
import networker.geo_math as gm
from networker.classes.geograph import GeoGraph
import warnings
import os

"""
Package for reading/writing networkx based GeoGraphs
Note:  these wrap existing networkx functions for custom behavior
"""

def load_shp(shp_path):
    """ loads a shapefile into a networkx based GeoGraph object

    Args:
        shp_path:  string path to a line or point shapefile

    Returns:
        geograph:  GeoGraph 

    """

    # NOTE:  if shp_path is unicode io doesn't work for some reason
    shp_path = shp_path.encode('ascii', 'ignore')
    g = nx.read_shp(shp_path)
    coords = dict(enumerate(g.nodes()))

    driver = ogr.GetDriverByName('ESRI Shapefile')
    shp = driver.Open(shp_path)
    layer = shp.GetLayer()

    spatial_ref = layer.GetSpatialRef()
    proj4 = None
    if not spatial_ref:
        if gm.is_in_lon_lat(coords):
           proj4 = gm.PROJ4_LATLONG
        else:
            warnings.warn("Spatial Reference could not be set for {}".\
            format(shp_path))

    else:
        proj4 = spatial_ref.ExportToProj4()

    g = nx.convert_node_labels_to_integers(g)

    return GeoGraph(srs=proj4, coords=coords, data=g)


def write_shp(geograph, shp_dir):
    """ writes a shapefile from a networkx based GeoGraph object

    Args:
        geograph: GeoGraph object
        shp_dir:  string path to dir to write shape files

    """

    assert geograph.is_aligned() 

    # looks like networkx wants us to relabel nodes by their coords
    tup_map = {i: tuple(coords) for i, coords in geograph.coords.items()}
    
    # copy geograph to plain networkx graph 
    # (relabeling a GeoGraph doesn't seem to work)
    nx_coord_graph = nx.Graph(data=geograph)
    nx.relabel_nodes(nx_coord_graph, tup_map, copy=False)

    nx.write_shp(nx_coord_graph, shp_dir)

    if geograph.srs:
        # write srs info to prj file (nx seems to miss this)
        sr = osr.SpatialReference()
        sr.ImportFromProj4(geograph.srs)
        main_prj_filename = shp_dir + '.prj'
        edge_prj_filename = os.path.join(shp_dir, 'edges.prj')
        node_prj_filename = os.path.join(shp_dir, 'nodes.prj')
        def write_prj(prj_filename):
            out = open(prj_filename, 'w')
            out.write(sr.ExportToWkt())
            out.close()
        
        write_prj(main_prj_filename)
        write_prj(edge_prj_filename)
        write_prj(node_prj_filename)
