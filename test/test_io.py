# -*- coding: utf-8 -*-

import networkx as nx
from osgeo import ogr
import operator
import networker.io as nio

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
