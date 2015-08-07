# -*- coding: utf-8 -*-

import numpy as np
import networkx as nx

from copy import deepcopy
from rtree import Rtree

from networker.classes.kdtree import KDTree
from networker.classes.unionfind import UnionFind, PriorityQueue
from networker.classes.geograph import GeoGraph

from networker.geomath import ang_to_vec_coords, \
                                  spherical_distance, \
                                  euclidean_distance, \
                                  make_bounding_box, \
                                  line_subgraph_intersection, \
                                  square_distance


def mod_boruvka(G, subgraphs=None, rtree=None):

    """
    algorithm to calculate the minimum spanning forest of nodes in GeoGraph G
    with 'budget' based restrictions on edges.

    Uses a modified version of Boruvka's algorithm

    NOTE:  subgraphs is modified as a side-effect...may remove in future
        (useful for testing right now)

    Args:
        G:  GeoGraph of nodes to be connected if appropriate
            Nodes should have 'budget' attribute
        subgraphs:  UnionFind data structure representing existing network's
            connected components AND the 'fake' nodes projected onto it.  This
            is the basis for the agglomerative nearest neighbor approach in
            this algorithm.
        rtree:  RTree based index of existing network

    Returns:
        GeoGraph: representing minimum spanning forest of G subject to the
            budget based restrictions
    """

    # special case (already MST)
    if G.number_of_nodes() < 2:
        return G

    V = set(G.nodes())
    coords = np.row_stack(G.coords.values())
    projcoords = ang_to_vec_coords(coords) if G.is_geographic() else coords
    kdtree = KDTree(projcoords)

    # Handle "dead" components
    D = set()

    if subgraphs is None:
        if rtree is not None:
            raise ValueError('RTree passed without UnionFind')

        rtree = Rtree()
        # modified to handle queues, children, mv
        subgraphs = UnionFind()

    # <helper_functions>
    def candidate_components(C):
        """
        return the set of candidate nearest components for the connected
        component containing C.  Do not consider those in C's connected
        component or those that are 'dead'.
        """
        component_set = subgraphs.component_set(C)
        return list((V - set(component_set)) - D)

    def update_nn_component(C, candidates):
        """
        find the nearest neighbor pair for the connected component
        represented by c.  candidates represents the list of
        components from which to select.
        """

        (v, vm) = subgraphs.queues[C].top()

        # vm ∈ C {not a foreign nearest neighbor}
        # go through the queue until an edge is found between this node
        # and the set of candidates, updating the neighbors in the connected
        # components queue in the process.
        while vm not in candidates:
            subgraphs.queues[C].pop()
            um, _ = kdtree.query_subset(projcoords[v], candidates)
            dm = square_distance(projcoords[v], projcoords[um])
            subgraphs.push(subgraphs.queues[C], (v, um), dm)
            # Note:  v will always be a vertex in this connected component
            #        vm *may* be external
            (v, vm) = subgraphs.queues[C].top()

        return (v, vm)

    def is_fake(node):
        """
        Tests whether the node is a projection on the existing grid
        """
        return subgraphs.budget[node] == np.inf

    # Test whether the component is dead
    # i.e. can never connect to another node
    def is_dead(c, nn_dist):
        return not is_fake(c) and subgraphs.budget[c] < nn_dist

    # "true" distance between components
    def component_dist(c1, c2):
        dist = 0
        if G.is_geographic():
            dist = spherical_distance([coords[c1], coords[c2]])
        else:
            dist = euclidean_distance([coords[c1], coords[c2]])
        return dist

    # </helper_functions>

    # Initialize the connected components holding a single node
    # and push the nearest neighbor into its queue
    for v in V:
        vm, _ = kdtree.query_subset(projcoords[v], list(V - {v}))
        dm = square_distance(projcoords[v], projcoords[vm])
        subgraphs.add_component(v, budget=G.node[v]['budget'])
        subgraphs.push(subgraphs.queues[v], (v, vm), dm)

        # Add to dead set if warranted
        nn_dist = component_dist(v, vm)

        if is_dead(v, nn_dist):
            # here components are single nodes
            # so no need to worry about adding children to dead set
            if v not in D:
                D.add(v)

    Et = []  # Initialize MST edges to empty list
    last_state = None

    # MST is complete when no progress was made in the prior iteration
    while Et != last_state:

        # This is a candidate list of edges that might be added to the MST
        Ep = PriorityQueue()

        # ∀ C of G; where c <- connected component
        for C in subgraphs.connected_components(component_subset=V):

            candidates = candidate_components(C)

            # Skip if no valid candidates
            if not candidates:
                continue

            (v, vm) = update_nn_component(C, candidates)

            # Add to dead set if warranted
            nn_dist = component_dist(v, vm)

            if is_dead(C, nn_dist):
                # add all dead components to the dead set D
                # (note that fake nodes can never be dead)
                for c in subgraphs.component_set(C):
                    if c not in D and not is_fake(c):
                        D.add(c)

        # One more round to root out connections to dead components
        # found in above iteration.
        # Need to do this BEFORE pushing candidate edges into Ep.
        # Otherwise we might be testing only 'dead' candidates
        # and therefore mistakenly think we were done (since
        # no new edges would have been added)
        for C in subgraphs.connected_components(component_subset=V):

            candidates = candidate_components(C)

            # Skip if no valid candidates
            if not candidates:
                continue

            (v, vm) = update_nn_component(C, candidates)

            # Calculate nn_dist for comparison to mv later
            nn_dist = component_dist(v, vm)

            # Append the top priority edge from the subgraph to the candidate
            # edge set
            Ep.push((v, vm, nn_dist), nn_dist)

        last_state = deepcopy(Et)

        # Candidate Test
        # At this point we have all of our nearest neighbor component edge
        # candidates defined for this "round"
        #
        # Now test all candidate edges in Ep for cycles and satisfaction of
        # custom criteria
        while Ep._queue:
            (um, vm, dm) = Ep.pop()

            # if doesn't create cycle
            # and subgraphs have enough MV
            # and we're not connecting 2 fake nodes
            # then allow the connection
            if subgraphs[um] != subgraphs[vm] and \
               (subgraphs.budget[subgraphs[um]] >= dm or is_fake(um)) and \
               (subgraphs.budget[subgraphs[vm]] >= dm or is_fake(vm)) and \
               not (is_fake(um) and is_fake(vm)):

                # doesn't create cycles from line segment intersection
                invalid_edge, intersections = \
                    line_subgraph_intersection(subgraphs, rtree,
                                               coords[um], coords[vm])

                if not invalid_edge:
                    # edges should not intersect a subgraph more than once
                    assert(filter(lambda n: n > 1,
                                  intersections.values()) == [])

                    # merge the subgraphs
                    subgraphs.union(um, vm, dm)

                    # For all intersected subgraphs update the mv to that
                    # created by the edge intersecting them,
                    # TODO: This should be updated in not such a naive method
                    map(lambda (n, _): subgraphs.union(um, n, 0),
                        filter(lambda (n, i): i == 1 and
                        subgraphs[n] != subgraphs[um],
                        intersections.iteritems()))

                    # index the newly added edge
                    box = make_bounding_box(coords[um], coords[vm])

                    # Object is (u.label, v.label), (u.coord, v.coord)
                    rtree.insert(hash((um, vm)), box,
                                 obj=((um, vm), (coords[um], coords[vm])))
                    Et += [(um, vm, {'weight': dm})]

    # create new GeoGraph with results
    result = G.copy()
    result.coords = G.coords
    result.remove_edges_from(result.edges())
    result.add_edges_from(Et)
    return result
