# -*- coding: utf-8 -*-

import numpy as np
import networkx as nx

from copy import deepcopy
from rtree import Rtree

from networker.classes.unionfind import UnionFind, PriorityQueue
from networker.classes.geograph import GeoGraph

from networker.geomath import ang_to_vec_coords, \
                                  spherical_distance, \
                                  euclidean_distance, \
                                  make_bounding_box, \
                                  line_subgraph_intersection, \
                                  square_distance


def mod_kruskal(G, subgraphs=None, rtree=None):

    """
    algorithm to compute the euclidean minimum spanning forest of nodes in
    GeoGraph G with 'budget' based restrictions on edges

    Uses a modified version of Kruskal's algorithm

    NOTE:  subgraphs is modified as a side-effect...may remove in future
        (useful for testing right now)

    Args:
        G:  GeoGraph of nodes to be connected if appropriate
            Nodes should have 'budget' attribute

        subgraphs:  UnionFind data structure representing existing network's
            connected components AND the 'fake' nodes projected onto it.  This
            is the basis for the agglomerative nearest neighbor approach in
            this algorithm.

            NOTE:  ONLY the existing networks components are represented in
            the subgraphs argument.  The nodes in G will be added within
            this function

        rtree:  RTree based index of existing network

    Returns:
        GeoGraph: representing minimum spanning forest of G subject to the
            budget based restrictions
    """

    # special case (already MST)
    if G.number_of_nodes() < 2:
        return G

    def is_fake(node):
        """
        Tests whether the node is a projection on the existing grid,
        using its MV
        """
        return subgraphs.budget[node] == np.inf

    # handy to have coords array
    coords = np.row_stack(G.coords.values())

    if subgraphs is None:
        assert rtree is not None, \
            "subgraphs (disjoint set) required when rtree is passed"

        rtree = Rtree()

        # modified to handle queues, children, mv
        subgraphs = UnionFind()

    # add nodes and budgets from G to subgraphs as components
    for node in G.nodes():
        subgraphs.add_component(node, budget=G.node[node]['budget'])

    # get fully connected graph
    g = G.get_connected_weighted_graph()

    # edges in MSF
    Et = []
    # connect the shortest safe edge until all edges have been tested
    # at which point, we have made all possible connections
    for u, v, w in sorted(g.edges(data=True), key=lambda x: x[2]['weight']):
        # if doesn't create cycle
        # and subgraphs have enough MV
        # and we're not connecting 2 fake nodes
        # then allow the connection
        w = w['weight']
        if subgraphs[u] != subgraphs[v] and \
           (subgraphs.budget[subgraphs[u]] >= w or is_fake(u)) and \
           (subgraphs.budget[subgraphs[v]] >= w or is_fake(v)) and \
           not (is_fake(u) and is_fake(v)):

            # doesn't create cycles from line segment intersection
            invalid_edge, intersections = \
                line_subgraph_intersection(subgraphs, rtree,
                                           coords[u], coords[v])

            if not invalid_edge:
                # edges should not intersect a subgraph more than once
                assert(filter(lambda n: n > 1,
                       intersections.values()) == [])

                # merge the subgraphs
                subgraphs.union(u, v, w)

                # For all intersected subgraphs update the mv to that
                # created by the edge intersecting them
                map(lambda (n, _): subgraphs.union(u, n, 0),
                                   filter(lambda (n, i): i == 1 and
                                          subgraphs[n] != subgraphs[u],
                                          intersections.iteritems()))

                # index the newly added edge
                box = make_bounding_box(coords[u], coords[v])

                # Object is (u.label, v.label), (u.coord, v.coord)
                rtree.insert(hash((u, v)), box,
                             obj=((u, v), (coords[u], coords[v])))
                Et += [(u, v, {'weight': w})]

    # create new GeoGraph with results
    result = G.copy()
    result.coords = G.coords
    result.remove_edges_from(result.edges())
    result.add_edges_from(Et)
    return result
