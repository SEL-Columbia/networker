# -*- coding: utf-8 -*-

__author__ = 'Brandon Ogle'

import numpy as np
import networkx as nx

from copy import deepcopy
from rtree import Rtree

from networkbuild.KDTree import KDTree
from networkbuild.utils import UnionFind, PriorityQueue, hav_dist, cartesian_projection, make_bounding_box, line_subgraph_intersection

def sq_dist(a,b):
    """Calculates square distance to reduce performance overhead of square root"""
    return np.sum((a-b)**2)

def modBoruvka(T, subgraphs=None, rtree=None, spherical_coords=True):

    # special case (already MST)
    if T.number_of_nodes() < 2:
        return T

    V = set(T.nodes(data=False))
    coords = np.row_stack(nx.get_node_attributes(T, 'coords').values())
    projcoords = cartesian_projection(coords) if spherical_coords else coords
    kdtree = KDTree(projcoords)

    # Handle "dead" components
    D = set()

    if subgraphs is None:
        if rtree != None: raise ValueError('RTree passed without UnionFind')

        rtree = Rtree()
        # modified to handle queues, children, mv
        subgraphs = UnionFind(T)

    # Tests whether the node is a projection on the existing grid, using its MV
    is_fake = lambda n: subgraphs.mv[n] == np.inf
    
    # Test whether the component is dead 
    # i.e. can never connect to another node
    def is_dead(c, nn_dist):
        return subgraphs.mv[c] < nn_dist

    #                ∀ v∈V
    # find the nearest neighbor for all nodes,
    # initialize a singleton subgraph and populate
    # its a queue with the nn edge where dist is priority
    for v in V:
        vm, _ = kdtree.query_subset(projcoords[v], list(V - {v}))
        dm = sq_dist(projcoords[v], projcoords[vm])
        C = subgraphs[v]
        subgraphs.push(subgraphs.queues[C], (v, vm), dm)

        # Add to dead set if warranted
        nn_dist = 0
        if spherical_coords:
            nn_dist = hav_dist(coords[v], coords[vm])
        else:
            nn_dist = np.sqrt(sq_dist(coords[v], coords[vm]))
        
        if is_dead(C, nn_dist):
            if C not in D: D.add(C)


    Et = [] # Initialize MST edges to empty list
    last_state = None

    # MST is complete when no progress was made in the prior iteration
    while Et != last_state:

        # This is an intermediary list of edges that might be added to the MST
        Ep = PriorityQueue()

        #∀ C of T; where C <- connected component
        for C in subgraphs.connected_components():

            (v, vm) = subgraphs.queues[C].top()
            component_set = subgraphs.component_set(v)
            djointVC = list(V - set(component_set))

            # if V ∈ C then minimum spanning tree, solved in last
            # iteration, continue to save state and terminate loop
            if not djointVC:
                continue

            # vm ∈ C {not a foreign nearest neighbor}
            # go through the queue until a edge is found that connects two subgraphs
            # while in the loop update the items in the queue,
            # preventing edges between nodes in the same subgraph
            while vm in component_set:
                subgraphs.queues[C].pop()
                um, _ = kdtree.query_subset(projcoords[v], djointVC)
                dm = sq_dist(projcoords[v], projcoords[um])
                subgraphs.push(subgraphs.queues[C], (v,um), dm)
                # Note:  v will always be a vertex in this connected component
                #        vm *may* be external
                (v,vm) = subgraphs.queues[C].top()

            # Add to dead set if warranted
            nn_dist = 0
            if spherical_coords:
                nn_dist = hav_dist(coords[v], coords[vm])
            else:
                nn_dist = np.sqrt(sq_dist(coords[v], coords[vm]))
            
            if is_dead(C, nn_dist):
                if C not in D:  D.add(C)

        # One more round to root out connections to dead components
        for C in subgraphs.connected_components():

            (v, vm) = subgraphs.queues[C].top()
            component_set = subgraphs.component_set(v)

            # Find nearest excluding Dead components too
            # Dead component sets
            dead_set = set()
            for c in D:
                dead_set = dead_set.union(set(subgraphs.component_set(c)))

            djointVC = list((V - set(component_set)) - dead_set)

            # if V ∈ C then minimum spanning tree, solved in last
            # iteration, continue to save state and terminate loop
            if not djointVC:
                continue

            # Look past both components of the sub-graph AND dead neighbors 
            while vm not in djointVC:
                subgraphs.queues[C].pop()
                um, _ = kdtree.query_subset(projcoords[v], djointVC)
                dm = sq_dist(projcoords[v], projcoords[um])
                subgraphs.push(subgraphs.queues[C], (v,um), dm)
                # v will always be one of the connected components
                # vm *may* be a a component external to these connected components
                (v,vm) = subgraphs.queues[C].top()

            # Calculate nn_dist for comparison to mv later 
            nn_dist = 0
            if spherical_coords:
                nn_dist = hav_dist(coords[v], coords[vm])
            else:
                nn_dist = np.sqrt(sq_dist(coords[v], coords[vm]))
 
            # Append the top priority edge from the subgraph to the intermediary edgelist
            Ep.push((v, vm, nn_dist), nn_dist)

        last_state = deepcopy(Et)

        # add all the edges in E' to Et so long as no cycles are created
        while Ep._queue:
            (um, vm, dm) = Ep.pop()
            # if doesn't create cycle and subgraph has enough MV
            if subgraphs[um] != subgraphs[vm] and (subgraphs.mv[subgraphs[um]] >= dm or is_fake(um)):
                # test that the connecting subgraph can receive the MV
                if subgraphs.mv[subgraphs[vm]] >= dm or is_fake(vm):

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
