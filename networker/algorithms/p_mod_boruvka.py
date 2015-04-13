__author__ = 'Brandon Ogle'

import numpy as np
import networkx as nx
import multiprocessing as mp

from copy import deepcopy
from rtree import Rtree

from networker.classes.kdtree import KDTree
from networker.classes.unionfind import UnionFind, PriorityQueue

from networker.geo_math import spherical_distance, \
                                  make_bounding_box, \
                                  line_subgraph_intersection, \
                                  ang_to_vec_coords

def p_mod_boruvka(G, subgraphs=None, rtree=None):
    
    V = set(T.nodes(data=False))
    coords = np.row_stack(nx.get_node_attributes(T, 'coords').values())
    projcoords = ang_to_vec_coords(coords)

    kdtree = KDTree(projcoords)

    if subgraphs is None:
        if rtree != None: raise ValueError('RTree passed without UnionFind')

        rtree = Rtree()
        # modified to handle queues, children, mv
        subgraphs = UnionFind()

    # Tests whether the node is a projection on the existing grid, using its MV
    is_fake = lambda n: subgraphs.budget[n] == np.inf

    def find_nn(node_tuple):
        u, up = node_tuple
        v, _ = kdtree.query_subset(up, list(V - {u}))
        return u, v, spherical_distance([coords[u], coords[v]])

    # find the nearest neigbor of all nodes
    p = mp.Pool(processes=6)
    neighbors = p.map(find_nn, enumerate(projcoords))
    p.close()

    # push the results into their respective queues
    for u, v, d in neighbors:
        subgraphs.add_component(u, budget=T.node[u]['budget'])
        subgraphs.queues[u].push((u, v), d)

    # list to hold mst edges
    Et = []
    last_state = None

    while Et != last_state:

        # consolidates top candidate edges from each subgraph
        Ep = PriorityQueue()

        def update_queues(component):

            q_top = subgraphs.queues[component].top()
            try:
                (u, v) = q_top
            except:
                return (None, None, None), np.inf

            component_set = subgraphs.component_set(u)
            disjointVC = list(V - set(component_set))

            if not disjointVC:
                return (None, None, None), np.inf

            while v in component_set:
                subgraphs.queues[component].pop()
                vprime, _ = kdtree.query_subset(projcoords[u], disjointVC)
                dm = spherical_distance([coords[u], coords[vprime]])
                subgraphs.queues[component].push((u, vprime), dm)
                (u, v) = subgraphs.queues[component].top()
            else:
                dm = spherical_distance([coords[u], coords[v]])

            return (u, v, dm), dm

        p = mp.Pool(processes=6)
        foreign_neighbors = map(update_queues, subgraphs.connected_components(component_subset=V))
        p.close()

        for neighbor in foreign_neighbors:
            obj, priority = neighbor
            if priority != np.inf:
                Ep.push(*neighbor)

        last_state = deepcopy(Et)
        # add all the edges in E' to Et so long as no cycles are created
        while Ep._queue:
            (um, vm, dm) = Ep.pop()
            # if doesn't create cycle and subgraph has enough MV
            if subgraphs[um] != subgraphs[vm] and (subgraphs.budget[subgraphs[um]] >= dm or is_fake(um)): 
                # test that the connecting subgraph can receive the MV
                if subgraphs.budget[subgraphs[vm]] >= dm or is_fake(vm):

                    # doesn't create cycles from line segment intersection
                    invalid_edge, intersections = line_subgraph_intersection(subgraphs, rtree, coords[um], coords[vm])

                    if not invalid_edge:
                        # valid edges should not intersect any subgraph more than once
                        assert(filter(lambda n: n > 1, intersections.values()) == [])

                        # merge the subgraphs
                        subgraphs.union(um, vm, dm)

                        # For all intersected subgraphs update the mv to that created by the edge intersecting them,
                        # TODO: This should be updated in not such a naive method
                        map(lambda (n, _): subgraphs.union(um, n, 0), 
                                filter(lambda (n, i): i == 1 and subgraphs[n] != subgraphs[um], intersections.iteritems()))

                        # index the newly added edge
                        box = make_bounding_box(coords[um], coords[vm])

                        # Object is in form of (u.label, v.label), (u.coord, v.coord)
                        rtree.insert(hash((um, vm)), box, obj=((um, vm), (coords[um], coords[vm])))
                        Et += [(um, vm)]  

    T.remove_edges_from(T.edges())
    T.add_edges_from(Et)
    return T
