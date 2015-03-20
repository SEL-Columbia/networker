# -*- coding: utf-8 -*-

__author__ = 'Brandon Ogle'

import numpy as np
import networkx as nx

from copy import deepcopy
from rtree import Rtree

from networkbuild.KDTree import KDTree
from networkbuild.utils import UnionFind, PriorityQueue, \
                               make_bounding_box,\
                               line_subgraph_intersection

from networkbuild.geo_math import spherical_projection, \
                                  spherical_distance_scalar

def sq_dist(a,b):
    """
    Calculates square distance to reduce performance overhead of square root
    """
    return np.sum((a-b)**2)


def modBoruvka(T, subgraphs=None, rtree=None, spherical_coords=True):

    # special case (already MST)
    if T.number_of_nodes() < 2:
        return T

    V = set(T.nodes(data=False))
    coords = np.row_stack(nx.get_node_attributes(T, 'coords').values())
    projcoords = spherical_projection(coords) if spherical_coords else coords
    kdtree = KDTree(projcoords)

    # Handle "dead" components
    D = set()

    if subgraphs is None:
        if rtree != None: raise ValueError('RTree passed without UnionFind')

        rtree = Rtree()
        # modified to handle queues, children, mv
        subgraphs = UnionFind(T)

    # <helper_functions> 
    def candidate_components(C):  
        """
        return the set of candidate nearest components for the connected component 
        containing c.  Do not consider those in C's connected component or 
        those that are 'dead'.
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
            dm = sq_dist(projcoords[v], projcoords[um])
            subgraphs.push(subgraphs.queues[C], (v,um), dm)
            # Note:  v will always be a vertex in this connected component
            #        vm *may* be external
            (v,vm) = subgraphs.queues[C].top()

        return (v, vm)

    # Tests whether the node is a projection on the existing grid, using its MV
    is_fake = lambda n: subgraphs.mv[n] == np.inf
    
    # Test whether the component is dead 
    # i.e. can never connect to another node
    def is_dead(c, nn_dist):
        return not is_fake(c) and subgraphs.mv[c] < nn_dist

    # "true" distance between components
    def component_dist(c1, c2):
        dist = 0
        if spherical_coords:
            dist = spherical_distance_scalar([coords[c1], coords[c2]])
        else:
            dist = np.sqrt(sq_dist(coords[c1], coords[c2]))
        return dist
 
    # </helper_functions>

    # Initialize the connected components holding a single node 
    # and push the nearest neighbor into its queue
    for v in V:
        vm, _ = kdtree.query_subset(projcoords[v], list(V - {v}))
        dm = sq_dist(projcoords[v], projcoords[vm])
        C = subgraphs[v]
        subgraphs.push(subgraphs.queues[C], (v, vm), dm)

        # Add to dead set if warranted
        nn_dist = component_dist(v, vm) 
       
        if is_dead(C, nn_dist):
            # here components are single nodes
            # so no need to worry about adding children to dead set
            if C not in D: D.add(C)

    Et = [] # Initialize MST edges to empty list
    last_state = None

    # MST is complete when no progress was made in the prior iteration
    while Et != last_state:

        # This is a candidate list of edges that might be added to the MST
        Ep = PriorityQueue()

        #∀ C of T; where c <- connected component
        for C in subgraphs.connected_components():

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
                    if c not in D and not is_fake(c):  D.add(c)

        # One more round to root out connections to dead components
        # found in above iteration.
        # Need to do this BEFORE pushing candidate edges into Ep.
        # Otherwise we might be testing only 'dead' candidates
        # and therefore mistakenly think we were done (since
        # no new edges would have been added)
        for C in subgraphs.connected_components():

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
        # Now test all candidate edges in Ep for cycles and satisfaction of custom
        # criteria
        while Ep._queue:
            (um, vm, dm) = Ep.pop()

            # if doesn't create cycle and subgraph has enough MV
            if subgraphs[um] != subgraphs[vm] and (subgraphs.mv[subgraphs[um]] >= dm or is_fake(um)):
                # test that the connecting subgraph can receive the MV AND
                # that both nodes are not fake
                # TODO:  Clean up this logic
                if (subgraphs.mv[subgraphs[vm]] >= dm or is_fake(vm)) and \
                   not (is_fake(um) and is_fake(vm)):

                    # doesn't create cycles from line segment intersection
                    invalid_edge, intersections = line_subgraph_intersection(subgraphs, rtree, coords[um], coords[vm])

                    if not invalid_edge:
                        # valid edges should not intersect any subgraph more than once
                        assert(filter(lambda n: n > 1, intersections.values()) == [])

                        # merge the subgraphs
                        subgraphs.union(um, vm, dm)

                        # For all intersected subgraphs update the mv to that created by
                        # the edge intersecting them,
                        # TODO: This should be updated in not such a naive method
                        map(lambda (n, _): subgraphs.union(um, n, 0),
                                filter(lambda (n, i): i == 1 and subgraphs[n] != subgraphs[um], intersections.iteritems()))

                        # index the newly added edge
                        box = make_bounding_box(coords[um], coords[vm])

                        # Object is in form of (u.label, v.label), (u.coord, v.coord)
                        rtree.insert(hash((um, vm)), box, obj=((um, vm), (coords[um], coords[vm])))
                        Et += [(um, vm, {'weight': dm})]


    T.remove_edges_from(T.edges())
    T.add_edges_from(Et)
    return T
