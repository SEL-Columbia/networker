# -*- coding: utf-8 -*-

import numpy as np
import networkx as nx
import networker.io as nio
from networker.classes.unionfind import UnionFind
import networker.geomath as gm

def test_project_point_2D():
    """
    Test whether point projection works in several important cases
    """

    # equivalence tolerance in distance squared
    SQ_TOLERANCE = 1e-10

    p1 = np.array([1.0, 1.0])
    l1 = np.array([[0.0, 0.0], [2.0, 0.0]])
    e = np.array([1.0, 0.0])

    p = gm.project_point_on_segment(p1, l1[0], l1[1])
    assert np.sum((e - p) ** 2) < SQ_TOLERANCE, \
        "projected point does not match expected"

    # now reverse the segment
    e = np.array([0.0, 0.0])
    p = gm.project_point_on_segment(p1, l1[0], -l1[1])
    assert np.sum((e - p) ** 2) < SQ_TOLERANCE, \
        "projected point does not match expected"

    # now with orthogonal segment
    e = np.array([0.0, 1.0])
    p = gm.project_point_on_segment(p1, l1[0], l1[1][[1, 0]])
    assert np.sum((e - p) ** 2) < SQ_TOLERANCE, \
        "projected point does not match expected"

    # now with orthogonal segment reversed
    e = np.array([0.0, 0.0])
    p = gm.project_point_on_segment(p1, l1[0], -l1[1][[1, 0]])
    assert np.sum((e - p) ** 2) < SQ_TOLERANCE, \
        "projected point does not match expected"


def test_segments_intersect():
    """
    Test various cases of intersection
    """

    def np_seg(seg):
        return map(np.array, seg)

    def check_intersects(intersect_fun):

        s1 = np_seg([[0.0, 0.0], [1.0, 0.0]])
        s2 = np_seg([[2.0, 0.0], [3.0, 0.0]])
        assert intersect_fun(s1[0], s1[1], s2[0], s2[1]) == False,\
            "segments {}, {} should not intersect".format(s1, s2)

        s1 = np_seg([[0.0, 0.0], [1.0, 1.0]])
        s2 = np_seg([[2.0, 0.0], [3.0, 0.0]])
        assert intersect_fun(s1[0], s1[1], s2[0], s2[1]) == False,\
            "segments {}, {} should not intersect".format(s1, s2)

        s1 = np_seg([[0.0, 0.0], [1.0, 0.0]])
        s2 = np_seg([[1.0, 0.0], [2.0, 0.0]])
        assert intersect_fun(s1[0], s1[1], s2[0], s2[1]) == True,\
            "segments {}, {} should intersect".format(s1, s2)

        s1 = np_seg([[0.0, 0.0], [1.0, 0.0]])
        s2 = np_seg([[0.5, 0.0], [2.0, 0.0]])
        assert intersect_fun(s1[0], s1[1], s2[0], s2[1]) == True,\
            "segments {}, {} should intersect".format(s1, s2)

        s1 = np_seg([[0.0, 0.0], [1.0, 0.0]])
        s2 = np_seg([[-1.0, 0.0], [2.0, 0.0]])
        assert intersect_fun(s1[0], s1[1], s2[0], s2[1]) == True,\
            "segments {}, {} should intersect".format(s1, s2)

    check_intersects(gm.segments_intersect)
    check_intersects(gm.segments_intersect_simple)


def test_project_point_3D():
    """
    Test whether 3D point projection works in several important cases
    """

    # equivalence tolerance in distance squared
    SQ_TOLERANCE = 1e-10

    # point 1/2 between x, y and z axes in northeastern hemisphere
    p1 = gm.ang_to_vec_coords([45.0, 45.0], radius=1.0)

    # v1, v2 represent arc between x and y axis
    v1 = np.array([1.0, 0.0, 0.0])
    v2 = np.array([0.0, 1.0, 0.0])

    # expect point to be 1/2 between x and y
    e = np.array([np.sin(np.pi/4), np.sin(np.pi/4), 0.0])

    p = gm.project_point_on_arc(p1, v1, v2, radius=1)
    assert np.sum((e - p) ** 2) < SQ_TOLERANCE, \
        "projected point on arc does not match expected"

    # point 1/2 between x, *negative* y and z axes in northeastern hemisphere
    p2 = gm.ang_to_vec_coords([-45.0, 45.0], radius=1.0)

    # expected point is v1 (i.e. x axis along unit sphere)
    e = v1

    p = gm.project_point_on_arc(p2, v1, v2, radius=1)
    assert np.sum((e - p) ** 2) < SQ_TOLERANCE, \
        "projected point not on arc does not match expected"


def test_spherical_dists():
    """
    compare dot vs haversine methods for spherical dist
    """

    xys = np.reshape(np.random.rand(20), (5,2,2))
    
    sum_diffs = np.abs(np.sum(gm.spherical_distance_dot(xys) - \
                gm.spherical_distance_haversine(xys)))

    assert sum_diffs < 1e-6, "dot and haversine results don't match"
    

def test_line_subgraph_intersection():
    """
    Test case where precision and odd geometry issues occur
    """
    # initialize network, nodes
    network = nio.load_shp("data/katsina/existing.shp", simplify=False)
    network.coords = {"g-" + str(n): network.coords[n] for n in network.nodes()}
    new_labels = ["g-" + str(n) for n in network.nodes()]
    nx.relabel_nodes(network, 
                     dict(zip(network.nodes(), new_labels)),
                     copy=False)
    nodes = nio.load_nodes("data/katsina/metrics.csv", "x", "y")
     
    # populate disjoint set of subgraphs
    subgraphs = UnionFind()
    # only one connected component, so just add all nodes associated
    # with first node
    net_nodes = network.nodes()
    parent = net_nodes[0]
    subgraphs.add_component(parent, budget=0)
    for node in net_nodes[1:]:
        subgraphs.add_component(node, budget=0)
        subgraphs.union(parent, node, 0)

    # now find projections onto grid
    rtree = network.get_rtree_index()
    projected = network.project_onto(nodes, rtree_index=rtree)
    projected.remove_nodes_from(network)

    assert len(projected.edges()) == 1, "should only be 1 projected edge"

    edge = projected.edges()[0]
    p1, p2 = projected.coords[edge[0]], projected.coords[edge[1]]

    invalid, subgraphs = gm.line_subgraph_intersection(subgraphs,
                                                    rtree,
                                                    p1, p2)

    assert not invalid, "edge should intersect network only once"
