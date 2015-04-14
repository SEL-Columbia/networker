# -*- coding: utf-8 -*-

import numpy as np
from networker.geomath import project_point_on_segment, \
                                  project_point_on_arc, \
                                  ang_to_vec_coords


def test_project_point_2D():
    """
    Test whether point projection works in several important cases
    """

    # equivalence tolerance in distance squared
    SQ_TOLERANCE = 1e-10

    p1 = np.array([1.0, 1.0])
    l1 = np.array([[0.0, 0.0], [2.0, 0.0]])
    e = np.array([1.0, 0.0])

    p = project_point_on_segment(p1, l1[0], l1[1])
    assert np.sum((e - p) ** 2) < SQ_TOLERANCE, \
        "projected point does not match expected"

    # now reverse the segment
    e = np.array([0.0, 0.0])
    p = project_point_on_segment(p1, l1[0], -l1[1])
    assert np.sum((e - p) ** 2) < SQ_TOLERANCE, \
        "projected point does not match expected"

    # now with orthogonal segment
    e = np.array([0.0, 1.0])
    p = project_point_on_segment(p1, l1[0], l1[1][[1, 0]])
    assert np.sum((e - p) ** 2) < SQ_TOLERANCE, \
        "projected point does not match expected"

    # now with orthogonal segment reversed
    e = np.array([0.0, 0.0])
    p = project_point_on_segment(p1, l1[0], -l1[1][[1, 0]])
    assert np.sum((e - p) ** 2) < SQ_TOLERANCE, \
        "projected point does not match expected"


def test_project_point_3D():
    """
    Test whether 3D point projection works in several important cases
    """

    # equivalence tolerance in distance squared
    SQ_TOLERANCE = 1e-10

    # point 1/2 between x, y and z axes in northeastern hemisphere
    p1 = ang_to_vec_coords([45.0, 45.0], radius=1.0)

    # v1, v2 represent arc between x and y axis
    v1 = np.array([1.0, 0.0, 0.0])
    v2 = np.array([0.0, 1.0, 0.0])

    # expect point to be 1/2 between x and y
    e = np.array([np.sin(np.pi/4), np.sin(np.pi/4), 0.0])

    p = project_point_on_arc(p1, v1, v2, radius=1)
    assert np.sum((e - p) ** 2) < SQ_TOLERANCE, \
        "projected point on arc does not match expected"

    # point 1/2 between x, *negative* y and z axes in northeastern hemisphere
    p2 = ang_to_vec_coords([-45.0, 45.0], radius=1.0)

    # expected point is v1 (i.e. x axis along unit sphere)
    e = v1

    p = project_point_on_arc(p2, v1, v2, radius=1)
    assert np.sum((e - p) ** 2) < SQ_TOLERANCE, \
        "projected point not on arc does not match expected"
