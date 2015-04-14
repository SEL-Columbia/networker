import numpy as np
import networkx as nx
from copy import deepcopy
from nose.tools import eq_

from networker.classes.unionfind import UnionFind
from networker.geo_math import spherical_distance

def init_network(n):

    coords = np.random.uniform(-10, 50, [n, 2])
    mv = np.random.uniform(7000000, 8000000, coords.shape[0])
    
    nodes = range(mv.size)

    
    graph = nx.Graph()
    graph.add_nodes_from(nodes)
    
    #add attributes
    nx.set_node_attributes(graph, 'coords', dict(enumerate(coords)))
    nx.set_node_attributes(graph,   'budget',   dict(enumerate(mv)))
    return graph
    
def test_union_budget():

    net = init_network(5000)
    subgraphs = UnionFind()
    nodes = net.nodes(data=True)
    pairs = zip(nodes[:-1], nodes[1:])

    mv = None
    for ((n1, d1), (n2, d2)) in pairs:
        if mv is None:
            mv = d1['budget'] 
        mv = (mv + d2['budget']) - \
             spherical_distance([d1['coords'], d2['coords']])

        subgraphs.add_component(n1, budget=d1['budget'])
        subgraphs.add_component(n2, budget=d2['budget'])
        d = spherical_distance([d1['coords'], d2['coords']])
        subgraphs.union(n1, n2, d)

    eq_(np.allclose(subgraphs.budget[subgraphs[1]], mv), True)


def test_component_functions():
    """
    Tests whether UnionFind component/connected_component methods
    work as expected

    """

    # demand nodes are 0-3 plus a 'fake node' at 4
    demand_components = set(range(5))
    # existing grid nodes
    external_components = set(['grid-0', 'grid-1', 'grid-2'])

    subgraphs = UnionFind()
    for i in demand_components:  subgraphs.add_component(i)
    for g in external_components:  subgraphs.add_component(g)

    # assign weights to eventual connected component roots
    subgraphs.weights[4] = np.inf
    subgraphs.weights[2] = np.inf
    subgraphs.weights['grid-1'] = np.inf

    # first connect the grid nodes and the fake node
    grid_union_pairs = [(4, 'grid-0'), ('grid-1', 'grid-2')]
    for g1, g2 in grid_union_pairs:
        subgraphs.union(g1, g2, 1)

    # should be 3 components at this point (2 within demand set)
    eq_(subgraphs.connected_components(component_subset=demand_components), set(range(5)))
    eq_(subgraphs.connected_components(), set(range(5) + ['grid-1']))

    # connect others (including a connection to the grid via fake node 4)
    union_pairs = [(0, 1), (2, 3), (0, 4)]
    for u1, u2 in union_pairs:
        subgraphs.union(u1, u2, 1)

    # test component sets
    eq_(set(subgraphs.component_set(0)), set([0, 1, 4, 'grid-0']))
    eq_(set(subgraphs.component_set(2)), set([2, 3]))

    # test connected components
    eq_(subgraphs.connected_components(), set([4, 2, 'grid-1']))
    # connected component with ('grid-1', 'grid-2') should be filtered out 
    eq_(subgraphs.connected_components(component_subset=demand_components), set([4, 2]))
     


