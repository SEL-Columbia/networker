# -*- coding: utf8 -*-
__author__ = 'Brandon Ogle'

import heapq
import osr
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
        self.neighborhood = {}
    
    def __getitem__(self, object):
        """Find and return the name of the set containing the object."""
        # check for previously unknown object
        if object not in self.parents:
            self.parents[object] = object
            self.weights[object] = 1
            self.mv[object] = self.graph.node[object]['mv']
            self.children[object] = [object]
            self.queues[object] = PriorityQueue()
            self.neighborhood[object] = set()

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
    

    def push(self, queue, item, priority):
        """Pushes an item into component queue, and updates the neighborhood"""
        u, v = item
        self.neighborhood[self[u]] |= {self[v]}
        self.neighborhood[self[v]] |= {self[u]}
        queue.push(item, priority)

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
        graphs = [g1, g2]; map(self.__getitem__, graphs)
        point_mv = np.array(map(self.mv.get, graphs))
        fake = point_mv == np.inf

        if all(fake):
            raise Exception('Path between fakes nodes')

        if any(fake):
            real = graphs[np.array([0, 1])[~fake]]
            grid = graphs[np.array([0, 1])[fake]]

            heaviest = self[real]
            smallest = self[grid]
        else:
            max_w = np.argmax(map(self.weights.get, graphs))
            heaviest = self[graphs[max_w]]
            smallest = self[graphs[~max_w]]

        self.weights[heaviest] += self.weights[smallest]
        self.parents[smallest] = heaviest

        self.children[heaviest] += self.children[smallest]
        self.children[smallest] = self.children[heaviest]

        self.queues[heaviest].merge(self.queues[smallest])
        self.queues[smallest] = self.queues[heaviest]

        if any(fake):
            if grid == smallest:
                self.mv[heaviest] -= d
                return

        self.mv[heaviest] += self.mv[smallest] - d
        self.mv[smallest] = self.mv[heaviest]

    def connected_components(self):
        """Return the roots for all disjoint sets"""
        return set([self.parents[r] for r in self.parents.keys() if not
                all('grid' in str(c) for c in self.children[self[r]])])

    def component_set(self, component):
        """Return the component set of the objects

        Args:
            n (obj): member node in one of the disjoint sets

        Returns:
            List of nodes in the same set as n
        """

        return [c for c in self.children[self[component]]
                if 'grid-' not in str(c)]

class PriorityQueue(object):

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
    Find the distance between a vector of (lat,lon) and
    the reference point (pcode_lat,pcode_lon).
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
    earth_radius = 6371010 # meters
    return haversine_angle * earth_radius


@jit
def cartesian_projection(coords):
    """projects x, y (lon, lat) coordinate pairs to 3D cartesian space"""
    
    rad = 6378137

    lon, lat = np.transpose(coords)
    cos_lat = np.cos(lat * np.pi / 180.0)
    sin_lat = np.sin(lat * np.pi / 180.0)
    cos_lon = np.cos(lon * np.pi / 180.0)
    sin_lon = np.sin(lon * np.pi / 180.0)
    xdim = rad * cos_lat * cos_lon
    ydim = rad * cos_lat * sin_lon
    zdim = rad * sin_lat
    return np.column_stack((xdim, ydim, zdim))

@jit
def make_bounding_box(u, v):
    stack = np.row_stack((u, v))
    xmin, xmax = np.argsort(stack[:, 0])
    ymin, ymax = np.argsort(stack[:, 1])
    return stack[xmin][0], stack[ymin][1], stack[xmax][0], stack[ymax][1]

@jit
def line_subgraph_intersection(subgraphs, rtree, p1, p2):
    """
    test for line segment intersection
    http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect

    TODO: Should compare this methodology to the dot product
    method to see which is more performant
    """

    box = make_bounding_box(p1, p2)

    # query for overlapping rectangles
    intersecting_bounds = rtree.intersection(box, objects=True)
    intersecting_subnets = defaultdict(int)

    # go through the possible intersections to validate
    for possible in intersecting_bounds:

        # Query object is in form of (u.label, v.label), (u.coord, v.coord)
        (up, vp), (p3, p4) = possible.object

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


    # Todo: If this edge is valid, we need to update
    # the mv for all intersecting subnets
    return False, intersecting_subnets

def project_point_to_segment(p, v1, v2):

    # Read line segment left to right
    if v1[0] > v2[0]:
        v1, v2 = v2, v1

    # Line segment vector (x, y)
    e1 = np.array([v2[0] - v1[0], v2[1] - v1[1]])
    # v1 -> p vector (x, y)
    e2 = np.array([p[0] - v1[0], p[1] - v1[1]])

    # Dot product and magnitude^2
    dp = np.dot(e1, e2)
    e_mag = np.sum(e1 ** 2)

    #P' of P on Line(v1, v2)
    proj = np.array([v1[0] + (dp * e1[0]) / e_mag,
                     v1[1] + (dp * e1[1]) / e_mag])

    # Constrain point at end points
    if proj[0] < v1[0]:
        return v1
    if proj[0] > v2[0]:
        return v2
    return proj

def csv_projection(path):
    """
    Args:
        path (str): path to the csv

    Returns:
        Projection (dict) of the csv coordinates if included in header
        else None
    """

    with open(path) as raw_text:
        header = raw_text.readline()

    if 'PROJ.4' in header:
        projection = string_to_proj4(header)


        return projection

def string_to_proj4(string):
    """
    Converts PROJ4 string to dict representation

    Args:
        string (str): PROJ4 string to convert

    Returns:
        proj4 (dict)

    """

    projection = {'proj': None,
                  'zone': None,
                  'ellps': None,
                  'datum': None,
                  'units': None}

    for key in projection.keys():
        start = string.find(key)
        if start != -1:
            start += len(key) + 1
            end = start + string[start:].find(' +')
            projection[key] = string[start:end]

    return projection


def sq_dist(a, b):

    return np.sum((a - b)**2, axis=1)


def all_dists(coords, spherical=True):
    """
    returns dist to nn and nn index as arrays indexed by coords index
    mainly used for testing
    """

    # get all perm's of coords
    a = np.tile(coords, (len(coords), 1))
    b = np.repeat(coords, len(coords), axis=0)
    all_dists = np.sqrt(sq_dist(a, b))
    if spherical:
        all_dists = get_hav_distance(a[:, 0], a[:, 1], b[:, 0], b[:, 1])

    zero_indices = np.array(range(len(coords))) * (len(coords) + 1)

    # so that mins are not the zero [i, i] vals
    all_dists[zero_indices] = np.inf
    full_dist_matrix = all_dists.reshape(len(coords), len(coords))
    return full_dist_matrix

def nn_dists(coords, spherical=True):

    full_dist_matrix = all_dists(coords, spherical)

    # find all minimum distances
    # apply min over ranges of the dist array
    min_dists = np.min(full_dist_matrix, axis=1)
    min_ind = np.argmin(full_dist_matrix, axis=1)

    return min_dists, min_ind


def utm_to_wgs84(coords, zone):
    """
    Converts an array of utm coords to latlon

    Args:
        coords (np.array): [[easting, northing]] coords to be converted
        zone   (int): utm zone associated with utm coords
    Returns:
        coords (np.array): [[Lon, Lat]]
    """
    utm_coordinate_system = osr.SpatialReference()

    # Set geographic coordinate system to handle lat/lon
    utm_coordinate_system.SetWellKnownGeogCS("WGS84")

    northing = coords[0, 1]
    is_northern = int(northing > 0)
    utm_coordinate_system.SetUTM(zone, is_northern)

    # Clone ONLY the geographic coordinate system
    lonlat_coordinate_system = utm_coordinate_system.CloneGeogCS()

    # create transformation helper
    converter = osr.CoordinateTransformation(utm_coordinate_system,
                                             lonlat_coordinate_system)

    return np.asarray(converter.TransformPoints(coords))[:, :-1]
