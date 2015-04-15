# -*- coding: utf8 -*-

import networkx as nx
import networker.geomath as gm 
import pyproj as prj
from networkx.readwrite import json_graph

# Utility functions
def csv_projection(path):
    """
    Get projection from csv file of geo data

    Args:
        path (str): path to the csv

    Returns:
        Proj4 projection if included in header else None
    """

    with open(path) as raw_text:
        header = raw_text.readline()

    if 'PROJ.4' in header:
        return header


# plot maps
def draw_np_graph(g, node_color='r', edge_color='b', node_label_field='ix',
                    edge_label_field=None, node_size=300):

    # transform to projected if not done so
    flat_coords = g.transform_coords(gm.PROJ4_FLAT_EARTH)

    node_pos = {nd: flat_coords[nd] for nd in g.nodes()}

    # handle main graph rendering
    if node_label_field:
        if node_label_field != 'ix':
            label_vals = nx.get_node_attributes(g, node_label_field)
            nx.draw_networkx(g, pos=node_pos, labels=label_vals,
                            node_color=node_color, edge_color=edge_color,
                            node_size=node_size)

        else: # label via ids
            nx.draw_networkx(g, pos=node_pos, with_labels=True,
                            node_color=node_color, edge_color=edge_color,
                            node_size=node_size)

    else:
        nx.draw_networkx(g, pos=node_pos, with_labels=False,
                        node_color=node_color, edge_color=edge_color,
                        node_size=node_size)

    # handle edge labels if needed
    if edge_label_field:
        edge_labels = nx.get_edge_attributes(g, edge_label_field)
        nx.draw_networkx_edge_labels(g, pos=node_pos, edge_labels=edge_labels)


# GeoGraph to js
def geograph_to_json(g):

    # transform to projected if not done so
    flat_coords = g.transform_coords(gm.PROJ4_FLAT_EARTH)

    g2 = g.copy()
    for nd in g.nodes():
        g2.node[nd]['coords'] = flat_coords[nd]

    js_g = json_graph.node_link_data(g2)
    return js_g
