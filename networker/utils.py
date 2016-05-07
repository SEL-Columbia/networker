# -*- coding: utf8 -*-

import networker.geomath as gm
import numpy as np
import networkx as nx
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


# GeoGraph to js
def geograph_to_json(g):

    # transform to projected if not done so
    flat_coords = g.transform_coords(gm.PROJ4_FLAT_EARTH)

    g2 = g.copy()
    for nd in g.nodes():
        g2.node[nd]['coords'] = flat_coords[nd]

    js_g = json_graph.node_link_data(g2)
    return js_g


def coords_dict_to_array2d(coords):
    """
    coords:  dict of coordinates

    returns 
    - numpy array of [num_coordsxcoord_dimension] 
    - index_node_id_map:  array of index->original_node_id
    """
    # create matrix placeholder
    num_rows = len(coords)
    num_cols = len(coords.values()[0])

    result_array = np.empty((num_rows, num_cols))
    index_node_id_map = []
    index = 0
    for key in coords.keys():
        index_node_id_map.append(key)
        result_array[index] = coords[key]
        index = index + 1
    
    return result_array, index_node_id_map


def array2d_to_coords_dict(array2d, index_node_id_map):
    """
    - array2d:  numpy array of [num_coordsxcoord_dimension] 
    - index_node_id_map:  array of index->original_node_id

    returns 
    coords:  dict of coordinate tuples with their original_node_id's
    """
    coords = dict()
    for i in range(array2d.shape[0]):
        key = index_node_id_map[i]
        coords[key] = tuple(array2d[i])

    return coords


def get_rounded_edge_sets(geograph, round_precision=8):
    """
    Get set of edges with coordinates rounded to a certain precision
    """
    tup_map = {i: 
               tuple(map(lambda c: round(c, round_precision), coords))
               for i, coords in geograph.coords.items()}

    edge_sets = map(lambda e: frozenset([tup_map[e[0]], tup_map[e[1]]]),
                    geograph.edges())
    return set(edge_sets)


def draw_geograph(g, node_color='r', edge_color='b', node_label_field=None,
                    edge_label_field=None, node_size=200):
    """
    Simple function to draw a geograph via matplotlib/networkx
    Uses geograph coords (projects if needed) as node positions
    """

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

"""
# spherical drawing helpers
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# globals for testing
global ax
def setup_3d_plot():
    global ax
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

def draw_wireframe(color="r"):

    global ax
    u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:50j]
    x=np.cos(u)*np.sin(v)
    y=np.sin(u)*np.sin(v)
    z=np.cos(v)
    ax.plot_wireframe(x, y, z, color=color)

def draw_arc(v1, v2, color="b", points_per_radian=100):
    # draw arc in 3d plot

    global ax
    arc_points = gm.get_arc_3D(v1, v2, points_per_radian=points_per_radian)
    ax.plot(arc_points[:,0], arc_points[:,1], arc_points[:,2], color=color)
"""
