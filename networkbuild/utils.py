# -*- coding: utf8 -*-
__author__ = 'Brandon Ogle'

import heapq
import numpy as np

from collections import defaultdict
from numba import jit

class Dict(dict):
    """dictionary allowing weakref"""
    pass

class UnionFind:

    def __init__(self, G): 
        """
        Creates and Empty UnionFind Structure:

        Args:
            G (NetworkX Graph)

        Notes:
            Union-find data structure also known as a Disjoint-set data structure

            Each unionFind instance X maintains a family of disjoint sets of
            hashable objects.

            This is a modified version of the data structure presented in Networkx
            nx.utils, with additonal properites to support the modified Boruvkas 
            algorithm. The original citations are listed below.

        NetworkX Citations
        ==================
        Union-find data structure. Based on Josiah Carlson's code,
        http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/215912
        with significant additional changes by D. Eppstein.
        http://www.ics.uci.edu/~eppstein/PADS/UnionFind.py

        """
        self.graph = G 
        self.weights = {}
        self.mv = {}
        self.parents = {}
        self.children = Dict() #This was previously used such that modifying it changed all refs pointing here
        self.queues = {}
    
    def __getitem__(self, object):
        """Find and return the name of the set containing the object."""

        # check for previously unknown object
        if object not in self.parents:
            self.parents[object] = object
            self.weights[object] = 1
            self.mv[object] = self.graph.node[object]['mv']
            self.children[object] = [object]
            self.queues[object] = PriorityQueue()
        
            return object

        # find path of objects leading to the root
        path = [object]
        root = self.parents[object]
        while root != path[-1]:
            path.append(root)
            root = self.parents[root]

        # compress the path and return
        for ancestor in path:
            self.parents[ancestor] = root
        return root
        
    def __iter__(self):
        """Iterate through all items ever found or unioned by this structure."""
        return iter(self.parents)

    def union(self, g1, g2, d):
        """
        Find the sets containing the objects and merge them all.
        During merge; children, mv and queues are all aggregated
        into the dominate or 'heaviest' set.
        
        Args:
            g1 (obj): Key of a member in the disjoint sets
            g2 (obj): Key of a member in the disjoint sets
             d (Num): distance between g1 and g2 
        """
        roots = [self[x] for x in [g1, g2]]
        heaviest = max([(self.weights[r],r) for r in roots])[1]
        r = [r for r in roots if r != heaviest][0]
        fakes = filter(lambda n: self.mv.get(n) == np.inf, roots)

        if fakes:
            if len(fakes) == 1 and fakes[0] == heaviest:
                heaviest, r = r, heaviest
            elif len(fakes) == 2: 
                raise Exception('Path Between Fake Nodes!')

        self.parents[r] = heaviest
        self.weights[heaviest] += self.weights[r]

        self.children[heaviest] += self.children[r]
        self.children[r] = self.children[heaviest]
    
        self.queues[heaviest].merge(self.queues[r])
        self.queues[r] = self.queues[heaviest]

        if fakes is None:
            self.mv[heaviest] = (self.mv[g1] + self.mv[g2]) - d
            self.mv[r] = self.mv[heaviest]

        # Leave fake node mv as infinity
        # The subgraph mv is then decreased by the distance to the fake node
        # Should this really be done? 
        else:
            self.mv[heaviest] = self.mv[heaviest] - d
    
    
    def connected_components(self):
        """Return the roots for all disjoint sets"""
        return set([self[r] for r in self.parents.keys()])
    
    def component_set(self, n):
        """Return the component set of the objects
        
        Args:
            n (obj): member node in one of the disjoint sets
        
        Returns: 
            List of nodes in the same set as n
        """
        
        return self.children[self[n]]

class PriorityQueue:
    
    def __init__(self):
        """
        Queue implementing highest-priority-in first-out.
        
        Note:
        Priority is cost based, therefore smaller values are prioritized
        over larger values. 
        """
        self._queue = []
        self._index = 0

    def push(self, item, priority):
        """
        Push an item into the queue.
        
        Args:
            item     (obj): Item to be stored in the queue
            priority (Num): Priority in which item will be retrieved from the queue
        """
        heapq.heappush(self._queue, (priority, self._index, item))
        self._index += 1

    def pop(self):
        """
        Removes the highest priority item from the queue

        Returns:
            obj: item with highest priority
        """
        return heapq.heappop(self._queue)[-1]
    
    def merge(self, other):
        """
        Given another queue, consumes each item in it
        and pushes the item and its priority into its own queue

        Args:
            other (PriorityQueue): Queue to be merged
        """
        while other._queue:
            priority,i,item = heapq.heappop(other._queue)
            self.push(item, priority)
        
    def top(self):
        """
        Allows peek at top item in the queue without removing it

        Returns:
            obj: if the queue is not empty otherwise None
        """
        try:
            return self._queue[0][-1]
        except:
            return None

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


@jit
def cartesian_projection(coords):
    """projects x, y (lon, lat) coordinate pairs to 3D cartesian space"""
    R = 6378137
    
    lon, lat = np.transpose(coords)
    cosLat = np.cos(lat * np.pi / 180.0)
    sinLat = np.sin(lat * np.pi / 180.0)
    cosLon = np.cos(lon * np.pi / 180.0)
    sinLon = np.sin(lon * np.pi / 180.0)
    rad = R
    x = rad * cosLat * cosLon
    y = rad * cosLat * sinLon
    z = rad * sinLat
    return np.column_stack((x,y,z))

@jit
def make_bounding_box(u, v):
    stack = np.row_stack((u, v))
    xmin, xmax = np.argsort(stack[:, 0])
    ymin, ymax = np.argsort(stack[:, 1])
    return stack[xmin][0], stack[ymin][1], stack[xmax][0], stack[ymax][1]

@jit
def line_subgraph_intersection(subgraphs, rtree, p1, p2):

    box = make_bounding_box(p1, p2)

    # query for overlapping rectangles
    intersecting_bounds = rtree.intersection(box, objects=True)
    intersecting_subnets = defaultdict(int)

    # go through the possible intersections to validate
    for possible in intersecting_bounds:
    
        # Query object is in form of (u.label, v.label), (u.coord, v.coord)
        (up, vp), (p3, p4) = possible.object            

        """
        test for line segment intersection
        http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect

        TODO: Should compare this methodology to the dot product method to see which is more performant
        """ 
        r = p2 - p1
        s = p4 - p3
        numerator = np.cross((p3 - p1), r)
        denominator = np.cross(r, s)
        
        if numerator == 0 and denominator == 0:
            # lines are colinear, test if overlapping 
            overlapping = (((p3[0]-p1[0]) < 0) != ((p3[0]-p2[0]) < 0)  != (
                            (p4[0]-p1[0]) < 0) != ((p4[0]-p2[0]) < 0)) or ((
                            (p3[1]-p1[1]) < 0) != ((p3[1]-p2[1]) < 0)  != (
                            (p4[1]-p1[1]) < 0) != ((p4[1] - p2[1]) < 0))
            if overlapping:
                # allow intersection if lines share an endpoint
                if (np.array_equal(p1, p3) and not np.array_equal(p2, p4)) or (
                    np.array_equal(p1, p4) and not np.array_equal(p2, p3)) or (
                    np.array_equal(p2, p3) and not np.array_equal(p1, p4)) or (
                    np.array_equal(p2, p4) and not np.array_equal(p1, p3)):
                    continue

                
                # Make sure something didn't go awry such that this edge doest represent a single subnet 
                assert(subgraphs[up] == subgraphs[vp])

                # Get the subgraph the segment intersects
                subgraph_parent = subgraphs[up]
                intersecting_subnets[subgraph_parent] += 1

                # If the subgraph is intersected in more than a single location 
                # this results in a 'cycle' and the segment is rejected
                if intersecting_subnets[subgraph_parent] > 1:
                    return True, intersecting_subnets

            else: continue
            
        if denominator == 0:
            # lines are parallel
            continue
            
        u = numerator / denominator
        t = np.cross((p3 - p1), s) / denominator
        
        intersecting = (0 <= t <= 1) and (0 <= u <= 1)
        if intersecting:
            # allow intersection if lines share an endpoint
            if (np.array_equal(p1, p4) and not np.array_equal(p2, p3)) or (
                np.array_equal(p2, p3) and not np.array_equal(p1, p4)) or (
                np.array_equal(p2, p4) and not np.array_equal(p1, p3)):
                continue
                
            # Make sure something went awry such that this edge doest represent a single subnet 
            assert(subgraphs[up] == subgraphs[vp])

            # Get the subgraph the segment intersects
            subgraph_parent = subgraphs[up]
            intersecting_subnets[subgraph_parent] += 1

            # If the subgraph is intersected in more than a single location 
            # this results in a 'cycle' and the segment is rejected
            if intersecting_subnets[subgraph_parent] > 1:
                return True, intersecting_subnets       
    

    # Todo: If this edge is valid, we need to update the mv for all intersecting subnets
    return False, intersecting_subnets

def project_point(L, q):
    """
    Linear system of equations to find projected point p of q on L
    http://cs.nyu.edu/~yap/classes/visual/03s/hw/h2/math.pdf
    """
    p0, p1 = L
    A = np.array([[p1[0] - p0[0], p1[1] - p0[1]],
                  [p0[1] - p1[1], p1[0] - p0[0]]])

    b = np.array([[-q[0]*(p1[0] - p0[0]) - q[1]*(p1[1] - p0[1])],
                  [-p0[1]*(p1[0] - p0[0]) + p0[0]*(p1[1]-p0[1])]])
    
    # A * [X,Y] = b
    # A^-1 * b = [X,Y]
    return np.ravel(np.linalg.solve(A, -b))
