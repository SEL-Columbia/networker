import networkx as nx

import numpy as np

from networker.algorithms.mod_boruvka import mod_boruvka
from networker.classes.unionfind import UnionFind
from networker.classes.geograph import GeoGraph
from rtree import Rtree


def graph_high_mvmax_long_edge():

    # This graph has an edge between components such that
    # if the min FNN is selected incorrectly, it will
    # produce a non minimal graph (this is just one case to test)
    # All vertices have sufficient MV

    # a(10) --- b(10) -------- c(10) --- d(10)
    #        3           7            3

    mv_max_values = [10, 10, 10, 10]
    coords = np.array([[0.0, 0.0], [3.0, 0.0], [10.0, 0.0], [13.0, 0.0]])
    coords_dict = dict(enumerate(coords))

    graph = GeoGraph(coords=coords_dict)
    nx.set_node_attributes(graph,   'budget',   dict(enumerate(mv_max_values)))

    return graph


def test_min_edges():

    # run mod_boruvka on graph with high mv max and long edge and make sure
    # that the result is an MST
    # of eachother (otherwise they should be connected)
    g = graph_high_mvmax_long_edge()

    subgraphs = UnionFind()
    rtree = Rtree()

    # build min span forest via mod_boruvka
    msf_g = mod_boruvka(g, subgraphs=subgraphs, rtree=rtree)

    # use networkx to build mst and compare
    coord_list = msf_g.coords.values()
    c = np.array(coord_list)
    all_dists = np.sqrt(((c[np.newaxis, :, :] - c[:, np.newaxis, :]) ** 2).
                    sum(2))

    complete_g = nx.Graph(all_dists)
    mst_g = nx.minimum_spanning_tree(complete_g)

    mst_edge_set = set([frozenset(e) for e in mst_g.edges()])
    msf_edge_set = set([frozenset(e) for e in msf_g.edges()])
    assert msf_edge_set == mst_edge_set
