import numpy as np
import networkx as nx
from copy import deepcopy
from nose.tools import eq_

from networkbuild.utils import UnionFind

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

def hav_dist(point1, point2):
    x1, y1 = point1
    x2, y2 = point2
    return get_hav_distance(y1, x1, y2, x2)

def get_hav_distance(lat, lon, pcode_lat, pcode_lon):
    """
    Find the distance between a vector of (lat,lon) and the reference point (pcode_lat,pcode_lon).
    """
    rad_factor = np.pi / 180.0  # degrees to radians for trig functions
    lat_in_rad = lat * rad_factor
    lon_in_rad = lon * rad_factor
    pcode_lat_in_rad = pcode_lat * rad_factor
    pcode_lon_in_rad = pcode_lon * rad_factor
    delta_lon = lon_in_rad - pcode_lon_in_rad
    delta_lat = lat_in_rad - pcode_lat_in_rad
    # Next two lines is the Haversine formula
    inverse_angle = (np.sin(delta_lat / 2) ** 2 + np.cos(pcode_lat_in_rad) *
                     np.cos(lat_in_rad) * np.sin(delta_lon / 2) ** 2)
    haversine_angle = 2 * np.arcsin(np.sqrt(inverse_angle))
    earth_radius =  6371010 # meters
    return haversine_angle * earth_radius

    
def TestUnionMV():

    net = init_network(5000)
    subgraphs = UnionFind(deepcopy(net))
    nodes = net.nodes(data=True)
    pairs = zip(nodes[:-1], nodes[1:])

    mv = None
    for ((n1, d1), (n2, d2)) in pairs:
        if mv is None:
            mv = d1['mv'] 
        mv = (mv + d2['mv']) - hav_dist(d1['coords'], d2['coords'])

        subgraphs[n1]; subgraphs[n2]
        d = hav_dist(d1['coords'], d2['coords'])
        subgraphs.union(n1, n2, d)

    eq_(np.allclose(subgraphs.mv[subgraphs[1]], mv), True)
