# -*- coding: utf-8 -*-

import ogr
import osr
import numpy as np
import pandas as pd
import networkx as nx
import networker.geomath as gm
from networker.utils import csv_projection
from networker.classes.geograph import GeoGraph
import warnings
import os

"""
Package for reading/writing networkx based GeoGraphs
Note:  these wrap existing networkx functions for custom behavior
"""


# TEMPORARY FIX:
# Copied from networkx library and modified, but I have a PR
# for this to be included in networkx lib
# https://github.com/networkx/networkx/pull/1491
def read_shp(path, simplify=True):
    """Generates a networkx.DiGraph from shapefiles. Point geometries are
    translated into nodes, lines into edges. Coordinate tuples are used as
    keys. Attributes are preserved, line geometries are simplified into start
    and end coordinates. Accepts a single shapefile or directory of many
    shapefiles.

    "The Esri Shapefile or simply a shapefile is a popular geospatial vector
    data format for geographic information systems software [1]_."

    Parameters
    ----------
    path : file or string
       File, directory, or filename to read.

    simplify:  Simplify line geometries to start and end coordinates
        If False, and line feature geometry has multiple segments, the
        attributes for that feature will be repeated for each edge comprising
        that feature

    Returns
    -------
    G : NetworkX graph

    Examples
    --------
    >>> G=nx.read_shp('test.shp') # doctest: +SKIP

    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/Shapefile
    """
    try:
        from osgeo import ogr
    except ImportError:
        raise ImportError("read_shp requires OGR: http://www.gdal.org/")

    if not isinstance(path, str):
        return

    net = nx.DiGraph()
    shp = ogr.Open(path)
    for lyr in shp:
        fields = [x.GetName() for x in lyr.schema]
        for f in lyr:
            flddata = [f.GetField(f.GetFieldIndex(x)) for x in fields]
            g = f.geometry()
            attributes = dict(zip(fields, flddata))
            attributes["ShpName"] = lyr.GetName()
            if g.GetGeometryType() == 1:  # point
                net.add_node((g.GetPoint_2D(0)), attributes)
            if g.GetGeometryType() == 2:  # linestring
                if simplify:
                    last = g.GetPointCount() - 1
                    attributes["Wkb"] = g.ExportToWkb()
                    attributes["Wkt"] = g.ExportToWkt()
                    attributes["Json"] = g.ExportToJson()
                    net.add_edge(g.GetPoint_2D(0), g.GetPoint_2D(last),
                                 attributes)
                else:
                    # NOTE: we do NOT add the geometry attributes (WKB...) here
                    # since they are not necessarily 1-1 with edges (and
                    # therefore can easily get out of sync on graph
                    # manipulation)
                    for i in range(0, g.GetPointCount() - 1):
                        net.add_edge(g.GetPoint_2D(i), g.GetPoint_2D(i+1),
                                     attributes)
    return net


def load_shp(shp_path, simplify=True):
    """ loads a shapefile into a networkx based GeoGraph object

    Args:
        shp_path:  string path to a line or point shapefile

        simplify:  Only retain start/end nodes of multi-segment lines

    Returns:
        geograph:  GeoGraph

    """

    # NOTE:  if shp_path is unicode io doesn't work for some reason
    shp_path = shp_path.encode('ascii', 'ignore')
    # TODO:  Call from networkx lib once PR has been
    # accepted and released
    # g = nx.read_shp(shp_path, simplify=simplify)
    g = read_shp(shp_path, simplify=simplify)
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
            warnings.warn("Spatial Reference could not be set for {}".
                          format(shp_path))

    else:
        proj4 = spatial_ref.ExportToProj4()

    g = nx.convert_node_labels_to_integers(g)

    return GeoGraph(srs=proj4, coords=coords, data=g)


def load_nodes(filename="metrics.csv", x_column="X", y_column="Y"):
    """
    load nodes csv into GeoGraph (nodes only)

    Args:
        filename:  nodal metrics csv file
        x_column, y_column:  col names to take x, y from

    Returns:
        GeoGraph of nodes including all attributes from input csv
    """

    input_proj = csv_projection(filename)
    header_row = 1 if input_proj else 0

    # read in the csv
    # NOTE:  Convert x,y via float cast to preserve precision of input
    metrics = pd.read_csv(filename, header=header_row,
                          converters={x_column: float, y_column: float})

    coord_cols = [x_column, y_column]
    assert all([hasattr(metrics, col) for col in coord_cols]), \
        "metrics file does not contain coordinate columns {}, {}".\
        format(x_column, y_column)

    # Stack the coords
    coords = np.column_stack(map(metrics.get, coord_cols))

    # set default projection
    if not input_proj:
        if gm.is_in_lon_lat(coords):
            input_proj = gm.PROJ4_LATLONG
        else:
            input_proj = gm.PROJ4_FLAT_EARTH

    coords_dict = dict(enumerate(coords))

    geo_nodes = GeoGraph(input_proj, coords_dict)

    # populate the rest of the attributes
    metrics_no_coords = metrics[metrics.columns.difference(coord_cols)]
    for row in metrics_no_coords.iterrows():
        index = row[0]
        attrs = row[1].to_dict()
        geo_nodes.node[index] = attrs

    return geo_nodes


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
