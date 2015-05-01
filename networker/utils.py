# -*- coding: utf8 -*-

import networkx as nx
import numpy as np
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
def draw_geograph(g, node_color='r', edge_color='b', node_label_field=None,
                    edge_label_field=None, node_size=200):

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
    """
    draw arc in 3d plot
    """

    global ax
    arc_points = gm.get_arc_3D(v1, v2, points_per_radian=points_per_radian)
    ax.plot(arc_points[:,0], arc_points[:,1], arc_points[:,2], color=color)
