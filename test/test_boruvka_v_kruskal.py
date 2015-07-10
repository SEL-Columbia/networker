# -*- coding: utf-8 -*-

import networkx as nx

import numpy as np

import networker.geomath as gm
from networker.algorithms.mod_boruvka import mod_boruvka
from networker.algorithms.mod_kruskal import mod_kruskal
from networker import networker_runner
from networker.classes.unionfind import UnionFind
from networker.classes.geograph import GeoGraph
from rtree import Rtree
import itertools


def nodes_plus_grid():
    """
    Return:  nodes as graph and grid as UnionFind/Rtree combo

    This example input demonstrates the "more" optimal nature 
    of mod_boruvka vs mod_kruskal.  

                          2(10)
                           |
            1(4)           |
  sqrt(5){ /|              | 
          / |              |
         /  | }3           | } 5
     (2)0   |              |
      1{|   |              |
    +-+-3-+-4-+-+-+-+-+-+-+5-+-+-+  <-- existing grid
  

    In this case, all nodes will be connected via either algorithm, 
    but the graph produced by mod_kruskal will have edge (2,5) whereas
    mod_boruvka will produce a graph with edge (0,1).  

    Therefore, the mod_boruvka graph is more optimal.  
    """

    mv_max_values = [2, 4, 10]
    coords = np.array([[0.0, 1.0], [1.0, 3.0], [10.0, 5.0]])
    coords_dict = dict(enumerate(coords))
    nodes = GeoGraph(gm.PROJ4_FLAT_EARTH, coords=coords_dict)

    nx.set_node_attributes(nodes, 'budget', dict(enumerate(mv_max_values)))

    grid_coords = np.array([[-5.0, 0.0], [15.0, 0.0]])
    grid = GeoGraph(gm.PROJ4_FLAT_EARTH, {'grid-' + str(n): c for n, c in
                    enumerate(grid_coords)})
    nx.set_node_attributes(grid, 'budget', {n: 0 for n in grid.nodes()})
    grid.add_edges_from([('grid-0', 'grid-1')])

    # now find projections onto grid
    rtree = grid.get_rtree_index()
    projected = grid.project_onto(nodes, rtree_index=rtree)
    projected.remove_nodes_from(grid)
    projected.remove_nodes_from(nodes)

    # populate disjoint set of subgraphs
    subgraphs = UnionFind()
    # only one connected component, so just associate all nodes 
    # with first node of grid
    parent = grid.nodes()[0]
    subgraphs.add_component(parent, budget=grid.node[parent]['budget'])
    for node in grid.nodes()[1:]:
        subgraphs.add_component(node, budget=grid.node[node]['budget'])
        subgraphs.union(parent, node, 0)

    # and the projected "fake" nodes
    for node in projected.nodes():
        subgraphs.add_component(node, budget=np.inf)
        subgraphs.union(parent, node, 0)

    # add projected nodes to node set
    nodes.add_nodes_from(projected, budget=np.inf)
    # merge coords
    nodes.coords = dict(nodes.coords, **projected.coords)

    return nodes, subgraphs, rtree 

# TODO:  Write test_optimality 
def test_optimality():

    # get test inputs, run mod_boruvka and mod_kruskal on them and 
    # ensure that mod_boruvka is more optimal

    # need copies of subgraphs, rtree because they are modified
    # by mod algos
    nodes1, subgraphs1, rtree1 = nodes_plus_grid()
    nodes2, subgraphs2, rtree2 = nodes_plus_grid()
    mod_b_result = mod_boruvka(nodes1, subgraphs1, rtree1) 
    mod_k_result = mod_kruskal(nodes2, subgraphs2, rtree2) 

    assert ((0, 1) in mod_b_result.edges() and
           (1, 4) not in mod_b_result.edges() and
           (0, 1) not in mod_k_result.edges() and
           (1, 4) in mod_k_result.edges()),\
           "mod_boruvka and mod_kruskal did not generate correct edges"

    def sum_edge_dists(geograph):
        tup_map = {i: tuple(coords) for i, coords in geograph.coords.items()}
        edge_coords = map(lambda e: [tup_map[e[0]], tup_map[e[1]]],
                        geograph.edges())
        dists = [gm.euclidean_distance(coord_pair)
                 for coord_pair in edge_coords]

        return sum(dists)

    assert sum_edge_dists(mod_b_result) < sum_edge_dists(mod_k_result),\
           "sum edge dists for mod_b should be < that for mod_k"
