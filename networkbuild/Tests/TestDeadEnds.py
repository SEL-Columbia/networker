import networkx as nx

import numpy as np

from networkbuild import modBoruvka
from networkbuild.utils import UnionFind
from rtree import Rtree
import itertools

def graph_with_dead_node():

    # a 'dead' node is one with insufficient mvMax to connect
    # to it's nearest neighbor (b in graph below)
    
    # a(10) ---- b(2) ----- c(10)
    #        4          4
    
    mv_max_values = [2, 10, 2, 10]
    coords = np.array([[-1, 0], [0, 0], [4, 0], [8, 0]])
  
    nodes = range(len(mv_max_values))
    graph = nx.Graph()
    graph.add_nodes_from(nodes)

    nx.set_node_attributes(graph, 'coords', dict(enumerate(coords)))
    nx.set_node_attributes(graph,   'mv',   dict(enumerate(mv_max_values)))

    return graph
  

def TestDeadBypass():

    # run modBoruvka on graph with dead node and make sure 
    # all connected components are NOT within their respective mvMax
    # of eachother (otherwise they should be connected)
    g = graph_with_dead_node()

    subgraphs = UnionFind(g)
    rtree = Rtree()

    # build min span forest via modBoruvka
    msf_g = modBoruvka(g, subgraphs=subgraphs, rtree=rtree, spherical_coords=False)

    # calculate all min distances between components
    # distances between all nodes (euclidean for now)
    coord_list = nx.get_node_attributes(msf_g, 'coords').values()
    c = np.array(coord_list)
    all_dists = np.sqrt(((c[np.newaxis, :, :] - c[:, np.newaxis, :]) ** 2).sum(2)) 
    
    # now get the component distances (min dist between components)
    components = subgraphs.connected_components()
    component_dists = {}

    # set init dists to inf
    for pair in itertools.product(components, components):
        component_dists[pair] = np.inf

    # keep the min distance between components
    for node_pair_dist in np.ndenumerate(all_dists):
        comp1 = subgraphs[node_pair_dist[0][0]]
        comp2 = subgraphs[node_pair_dist[0][1]]
        dist = node_pair_dist[1]
        if dist < component_dists[(comp1, comp2)]:
            component_dists[(comp1, comp2)] = dist
    
    # now check whether components are within 
    # their respective mvmax of eachother 
    # (if so, we've got a problem)
    missed_connections = []
    for pair in itertools.product(components, components):
        
        if pair[0] != pair[1] and \
           subgraphs.mv[pair[0]] >= component_dists[pair] and \
           subgraphs.mv[pair[1]] >= component_dists[pair]: 
            missed_connections.append(pair)

    assert len(missed_connections) == 0, "missed connections: " + str(missed_connections)








