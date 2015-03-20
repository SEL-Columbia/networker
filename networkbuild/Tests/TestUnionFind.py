import numpy as np
import networkx as nx
from copy import deepcopy
from nose.tools import eq_

from networkbuild.utils import UnionFind
from networkbuild.geo_math import spherical_distance_scalar

def init_network(n):

    coords = np.random.uniform(-10, 50, [n, 2])
    mv = np.random.uniform(7000000, 8000000, coords.shape[0])
    
    nodes = range(mv.size)

    
    graph = nx.Graph()
    graph.add_nodes_from(nodes)
    
    #add attributes
    nx.set_node_attributes(graph, 'coords', dict(enumerate(coords)))
    nx.set_node_attributes(graph,   'mv',   dict(enumerate(mv)))
    return graph
    
def TestUnionMV():

    net = init_network(5000)
    subgraphs = UnionFind(deepcopy(net))
    nodes = net.nodes(data=True)
    pairs = zip(nodes[:-1], nodes[1:])

    mv = None
    for ((n1, d1), (n2, d2)) in pairs:
        if mv is None:
            mv = d1['mv'] 
        mv = (mv + d2['mv']) - \
             spherical_distance_scalar([d1['coords'], d2['coords']])

        subgraphs[n1]; subgraphs[n2]
        d = spherical_distance_scalar([d1['coords'], d2['coords']])
        subgraphs.union(n1, n2, d)

    eq_(np.allclose(subgraphs.mv[subgraphs[1]], mv), True)
