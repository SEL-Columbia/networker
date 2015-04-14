# -*- coding: utf8 -*-

import networkx as nx


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

    from mpl_toolkits.basemap import Basemap
    m = Basemap(
            projection='merc',
            ellps='WGS84',
            llcrnrlon=0,
            llcrnrlat=0,
            urcrnrlon=1,
            urcrnrlat=1,
            lat_ts=0,
            resolution='i',
            suppress_ticks=True)

    node_pos = {nd: m(g.node[nd]['coords'][0], g.node[nd]['coords'][1])
                for nd in g.nodes()}

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
def network_to_json(g):

    from mpl_toolkits.basemap import Basemap
    from networkx.readwrite import json_graph
    m = Basemap(
            projection='merc',
            ellps='WGS84',
            llcrnrlon=0,
            llcrnrlat=0,
            urcrnrlon=1,
            urcrnrlat=1,
            lat_ts=0,
            resolution='i',
            suppress_ticks=True)

    g2 = g.copy()
    for nd in g.nodes():
        g2.node[nd]['coords'] = m(g.coords[nd][0], g.coords[nd][1])

    js_g = json_graph.node_link_data(g2)
    return js_g
