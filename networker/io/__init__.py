# -*- coding: utf-8 -*-

import ogr
import osr
import json
import numpy as np
import pandas as pd
import networkx as nx
import networker.geomath as gm
from networkx.readwrite import json_graph
from networker.utils import csv_projection
from networker.utils import nested_dict_getter
from networker.classes.geograph import GeoGraph
import warnings
import os

"""
Package for reading/writing networkx based GeoGraphs
Note:  these wrap existing networkx functions for custom behavior
"""


def read_shp(path, simplify=True, geom_attrs=True):
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

    simplify:  bool
        If ``True``, simplify line geometries to start and end coordinates.
        If ``False``, and line feature geometry has multiple segments, the
        non-geometric attributes for that feature will be repeated for each
        edge comprising that feature.

    geom_attrs: bool
        If ``True``, include the Wkb, Wkt and Json geometry attributes with
        each edge.

        NOTE:  if these attributes are available, write_shp will use them
        to write the geometry.  If nodes store the underlying coordinates for
        the edge geometry as well (as they do when they are read via
        this method) and they change, your geomety will be out of sync.


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
            # Note:  Using layer level geometry type
            if g.GetGeometryType() == ogr.wkbPoint:
                net.add_node((g.GetPoint_2D(0)), attributes)
            elif g.GetGeometryType() in (ogr.wkbLineString,
                                         ogr.wkbMultiLineString):
                for edge in edges_from_line(g, attributes, simplify,
                                            geom_attrs):
                    net.add_edge(*edge)
            else:
                raise ImportError("GeometryType {} not supported".
                                  format(g.GetGeometryType()))

    return net


def edges_from_line(geom, attrs, simplify=True, geom_attrs=True):
    """
    Generate edges for each line in geom
    Written as a helper for read_shp

    Parameters
    ----------

    geom:  ogr line geometry
        To be converted into an edge or edges

    attrs:  dict
        Attributes to be associated with all geoms

    simplify:  bool
        If ``True``, simplify the line as in read_shp

    geom_attrs:  bool
        If ``True``, add geom attributes to edge as in read_shp


    Yields:
        edges:  generator of edges
            each edge is a tuple of form
            (node1_coord, node2_coord, attribute_dict)
            suitable for expanding into a networkx Graph add_edge call
    """
    if geom.GetGeometryType() == ogr.wkbLineString:
        if simplify:
            edge_attrs = attrs.copy()
            last = geom.GetPointCount() - 1
            if geom_attrs:
                edge_attrs["Wkb"] = geom.ExportToWkb()
                edge_attrs["Wkt"] = geom.ExportToWkt()
                edge_attrs["Json"] = geom.ExportToJson()
            yield (geom.GetPoint_2D(0), geom.GetPoint_2D(last), edge_attrs)
        else:
            for i in range(0, geom.GetPointCount() - 1):
                pt1 = geom.GetPoint_2D(i)
                pt2 = geom.GetPoint_2D(i + 1)
                edge_attrs = attrs.copy()
                if geom_attrs:
                    segment = ogr.Geometry(ogr.wkbLineString)
                    segment.AddPoint_2D(pt1[0], pt1[1])
                    segment.AddPoint_2D(pt2[0], pt2[1])
                    edge_attrs["Wkb"] = segment.ExportToWkb()
                    edge_attrs["Wkt"] = segment.ExportToWkt()
                    edge_attrs["Json"] = segment.ExportToJson()

                yield (pt1, pt2, edge_attrs)

    elif geom.GetGeometryType() == ogr.wkbMultiLineString:
        for i in range(geom.GetGeometryCount()):
            geom_i = geom.GetGeometryRef(i)
            for edge in edges_from_line(geom_i, attrs, simplify, geom_attrs):
                yield edge


def load_shp(shp_path, simplify=True):
    """ 
    loads a shapefile into a networkx based GeoGraph object

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
    g = read_shp(shp_path, simplify=simplify, geom_attrs=False)
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
    """ 
    writes a shapefile from a networkx based GeoGraph object

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


def load_json(stream):
    """
    Args:
        stream: Open stream containing js

    Assumes the js is in networkx link-node format
    """
    
    js = json.load(stream)
    g = json_graph.node_link_graph(js)

    assert all([nd.has_key('coords') for nd in g.node.values()]),\
           "json node-link graph must have nodes with coords for GeoGraph"

    # get coords
    coords = [v['coords'] for v in g.node.values()]

    # set default projection
    input_proj = ""
    if gm.is_in_lon_lat(coords):
        input_proj = gm.PROJ4_LATLONG
    else:
        input_proj = gm.PROJ4_FLAT_EARTH

    coords_dict = {k: v['coords'] for k, v in g.node.items()}
    # now get rid of 'coords' key,val for each node
    for node in g.node.values():
        node.pop('coords', None)

    geo_nodes = GeoGraph(srs=input_proj, coords=coords_dict, data=g)
    return geo_nodes


def write_json(geograph, stream):
    """
    Args:
        geograph:  A GeoGraph object
        stream:  stream to output networkx link-node format json rep
    """
    g2 = geograph.copy()
    for nd in geograph.nodes():
        if isinstance(geograph.coords[nd], np.ndarray):
            g2.node[nd]['coords'] = geograph.coords[nd].tolist()
        else:
            g2.node[nd]['coords'] = geograph.coords[nd]

    js_g = json_graph.node_link_data(g2)
    stream.write(json.dumps(js_g))


def geojson_get_nodes(features):
    """
    Get all features that comply with GeoGraph node format within GeoJSON
    features collection

    Returns:
        nodes:  list of (node_id, node_attr_dict) tuples
        coords:  coordinates dict indexed by node_id
    """
    dict_getter = nested_dict_getter()
    nodes = []
    coords = {}
    for feature in features:
        # Only allow Point type as nodes for now
        if dict_getter(feature, ['geometry', 'type']) != 'Point':
            continue

        node_id = dict_getter(feature, ['properties', 'node_id'])
        # Nodes are only valid with a node_id property
        if node_id is None:
            continue

        node_attrs = feature['properties']
        node_attrs.pop('node_id')
        nodes.append((node_id, node_attrs))
        # At this point, we assume we have coordinates
        coords[node_id] = feature['geometry']['coordinates']

    return nodes, coords


def geojson_get_edges(features):
    """
    Get all features that comply with GeoGraph edge format within GeoJSON
    features collection

    Returns:
        edges:  list of (from_node_id, to_node_id, edge_attr_dict) tuples
    """
    dict_getter = nested_dict_getter()
    edges = []
    for feature in features:
        # Only allow LineString type as edges for now
        if dict_getter(feature, ['geometry', 'type']) != 'LineString':
            continue

        node_from_id = dict_getter(feature, ['properties', 'node_from'])
        node_to_id = dict_getter(feature, ['properties', 'node_to'])
        # Edges are only valid with node_from_id and node_to_id properties
        if node_from_id is None or node_to_id is None:
            continue

        edge_attrs = feature['properties']
        edge_attrs.pop('node_from')
        edge_attrs.pop('node_to')
        edges.append((node_from_id, node_to_id, edge_attrs))

    return edges

def load_geojson(geojson_path, directed=False):
    """
    Load GeoGraph from geojson

    See:  docs/geograph_geojson.md for format details

    Args:
        geojson_path:  Path to geojson

    Returns:
        geograph:  GeoGraph object
    """

    geojson = json.load(open(geojson_path, 'r'))
    dict_getter = nested_dict_getter()

    # build GeoGraph from a plain-old networkx graph
    g = nx.Graph()

    # Only add crs if it's explicitly set in crs.properties.name
    crs = dict_getter(geojson, ['crs', 'properties', 'name'])
    if crs:
        try:
            prj.Proj(crs)
        except Exception as e:
            raise SpatialReferenceInvalidException(
                "Spatial reference must comply with proj4, got {}" % crs)

    # TODO:  Apply geojson schema validation
    features = geojson['features']

    # Currently this separates nodes and coordinates due to structure of
    # GeoGraph.  This may change to allow Nodes and Edges to have distinct 
    # Geometry.
    nodes, coords = geojson_get_nodes(features)
    g.add_nodes_from(nodes)
    
    edges = geojson_get_edges(features)
    g.add_edges_from(edges)

    if crs is None and gm.is_in_lon_lat(coords.values()):
        crs = gm.PROJ4_LATLONG
    
    if crs is None:
        warnings.warn("Spatial Reference could not be set for {}".
                          format(geojson_path))

    # TODO:  More srs/projection validation?
    return GeoGraph(srs=crs, coords=coords, data=g)


def node_geojson():
    return  {
                "type": "Feature",
                "geometry": 
                {
                    "type": "Point",
                    "coordinates": []
                },
                "properties":
                {
                    "node_id": None
                }
            }

def edge_geojson():
    return  {
                "type": "Feature",
                "geometry": 
                {
                    "type": "LineString",
                    "coordinates": []
                },
                "properties":
                {
                    "node_from": None,
                    "node_to": None
                }
            }

def geojson_instance():
    return {
                "type": "FeatureCollection",
                "features": []
            }

def to_geojson(geograph):
    """
    Args:
        geograph:  A GeoGraph object

    Returns:
        geojson:  dict representing GeoGraph as dictionary
    """
    def node_to_geojson(node_id):
        node_geo_dict = node_geojson()
        node_geo_dict["properties"]["node_id"] = node_id
        node_geo_dict["properties"].update(geograph.node[node_id])
        node_geo_dict["geometry"]["coordinates"] = geograph.coords[node_id]
        return node_geo_dict

    def edge_to_geojson(node_from, node_to):
        edge_geo_dict = edge_geojson()
        edge_geo_dict["properties"]["node_from"] = node_from
        edge_geo_dict["properties"]["node_to"] = node_to
        edge_dict = geograph.edge[node_from][node_to]
        coordinates = [geograph.coords[node_from], geograph.coords[node_to]]
        # TODO:  is there a better way to handle custom coordinates than 
        # looking for "coordinates" attribute?
        if edge_dict.has_key("coordinates"):    
            coordinates = edge_dict.pop("coordinates")
        edge_geo_dict["properties"].update(edge_dict)
        edge_geo_dict["geometry"]["coordinates"] = coordinates
        return edge_geo_dict


    node_list = [node_to_geojson(node) for node in geograph.nodes()]
    edge_list = [edge_to_geojson(node_from, node_to) for 
                 node_from, node_to in geograph.edges()]
    geojson = geojson_instance()
    geojson["features"] = node_list + edge_list
    return geojson


def write_geojson(geograph, geojson_path):
    """
    Args:
        geograph:  A GeoGraph object
        geojson_path:  file to output GeoGraph GeoJSON format rep
    """
 
    geojson = to_geojson(geograph) 
    with open(geojson_path, 'w') as geojson_stream:
        geojson_stream.write(json.dumps(geojson))
