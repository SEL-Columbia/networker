import networkx as nx

from nose.tools import eq_

from networkbuild.classes import GeoGraph
from networkbuild.utils import UnionFind
import networkbuild.geo_math as gm
                                  
# fixtures
def network_nodes_projections():
    """
    Create a network and node graph to be merged, and
    the expected result of project_onto for testing
    """

    def to_dict(a):
        return {i: val for i, val in enumerate(a)}

    net_coords = [[0, 0], [3, 0], [3, 3]] 
    net_edges = [(0, 1), (1, 2)]
    node_coords = [[-2, 0], [1, 1], [4, -1], [4, 1]]

    projected_coords = [[0, 0], [1, 0], [3, 0], [3, 1]]
    
    g_net = GeoGraph(gm.PROJ4_FLAT_EARTH, to_dict(net_coords), data=net_edges)
    g_nodes = GeoGraph(gm.PROJ4_FLAT_EARTH, to_dict(node_coords))

    return g_net, g_nodes, projected_coords
    
def test_project_onto():

    net, nodes, projections = network_nodes_projections()

    new_net = net.project_onto(nodes)

    # check that the projected coords are in the new net
    assert all(p in new_net.coords.items() for p in projections), \
           "expected projection does not match actual"


