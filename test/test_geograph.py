# -*- coding: utf-8 -*-

import networkx as nx
import numpy as np

from networker.classes.geograph import GeoGraph
import networker.geomath as gm


# fixtures
def network_nodes_projections():
    """
    Create a network and node graph to be merged, and
    the expected result of project_onto for testing
    """

    net_coords = [[0.0, 0.0], [3.0, 0.0], [3.0, 3.0]]
    net_edges = [(0, 1), (1, 2)]
    node_coords = {3: [-2.0, 0.0], 4: [1.0, 1.0],
                   5: [4.0, -1.0], 6: [4.0, 1.0]}

    projected_coords = {3: [0.0, 0.0], 4: [1.0, 0.0],
                        5: [3.0, 0.0], 6: [3.0, 1.0]}

    g_net = GeoGraph(gm.PROJ4_FLAT_EARTH, dict(enumerate(net_coords)), data=net_edges)
    g_nodes = GeoGraph(gm.PROJ4_FLAT_EARTH, node_coords)

    return g_net, g_nodes, projected_coords


def test_project_onto():

    net, nodes, projections = network_nodes_projections()

    new_net = net.project_onto(nodes)

    # check that the projected coords are in the new net
    neighbor_coords = {node: [list(new_net.coords[neigh])
        for neigh in new_net.neighbors(node)]
        for node in projections.keys()}

    # need to 'cast' to list in case coords are numpy arrays
    tests = [projections[p] in neighbor_coords[p] for p in projections.keys()]
    assert all(tests), \
        "expected projection does not match actual"

    # ensure that the rtree based projection results in the same network
    rtree = net.get_rtree_index()
    new_net_rt = net.project_onto(nodes, rtree_index=rtree)

    assert nx.is_isomorphic(new_net, new_net_rt) and \
        [list(c) for c in new_net.coords.values()] == \
        [list(c) for c in new_net_rt.coords.values()], \
        "project_onto with and without rtree inputs don't match"


def test_connected_graph():
    g, _, _ = network_nodes_projections()

    g_conn = g.get_connected_weighted_graph()

    nodes = g.nodes()
    edges = set(zip(np.tile(nodes, len(nodes)), np.repeat(nodes, len(nodes))))
    diagonal = set(zip(nodes, nodes))
    edges = edges - diagonal 
    # order does NOT matter
    edges = set([frozenset(edge) for edge in edges])
    def dist_for_edge(edge):
        return gm.euclidean_distance([g.coords[edge[0]], g.coords[edge[1]]])

    weights = {edge: dist_for_edge(edge) for edge in \
                [tuple(set_edge) for set_edge in edges]}

    g_conn_weights = {(e[0],e[1]): e[2]['weight'] 
                        for e in g_conn.edges(data=True)}
    
    assert weights == g_conn_weights,\
        "fully connected edges/weights are not correct"
