# -*- coding: utf-8 -*-

import os
import networkx as nx
from osgeo import ogr
import operator
import networker.io as nio
import networker.geomath as gm
from networker.classes.geograph import GeoGraph

def test_edges_from_line():

    multiline = ogr.Geometry(ogr.wkbMultiLineString)
    lines = [
             ((0.0, 1.0), (1.0, 1.0)),
             ((2.0, 2.0), (3.0, 2.0))
            ]

    for line in lines:
        geom = ogr.Geometry(ogr.wkbLineString)
        for point in line:
            geom.AddPoint_2D(*point)

        multiline.AddGeometry(geom)

    attrs = {'name': 'MultiLine'}

    g = nx.Graph()
    for edge in nio.edges_from_line(multiline,
                                    attrs,
                                    simplify=False,
                                    geom_attrs=False):
        g.add_edge(*edge)

    g_expected = nx.Graph()
    g_expected.add_edges_from(lines, name='MultiLine')

    assert nx.is_isomorphic(g, g_expected,
                            node_match=operator.eq,
                            edge_match=operator.eq),\
           "expected graphs to be equal"

def test_load_write_json():
    """
    ensure that reading/writing js 'node-link' format works
    """

    os.mkdir(os.path.join('test', 'tmp'))
    node_dict = {0: [0,0], 1: [0,1], 2: [1,0], 3: [1,1]}
    g = GeoGraph(gm.PROJ4_LATLONG, node_dict)
    g.add_edges_from([(0,1),(1,2),(2,3)])
    json_file_path = os.path.join('test', 'tmp', 'g.js')
    nio.write_json(g, open(json_file_path, 'w'))

    g2 = nio.load_json(open(json_file_path, 'r'))
    os.remove(json_file_path)
    os.rmdir(os.path.join('test', 'tmp'))
    assert nx.is_isomorphic(g, g2,
                            node_match=operator.eq,
                            edge_match=operator.eq),\
           "expected written and read graphs to match"

def test_load_write_geojson():
    """
    ensure that reading geojson works
    """
    g = nio.load_geojson(os.path.join('data', 'geojson_sample.json'))
    print(g.edge)
    assert g.edge[0][1] == {'name': 'edge-0'}
    assert len(g.coords) == 2

