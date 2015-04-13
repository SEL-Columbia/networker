# -*- coding: utf-8 -*-

import osr
import networkx as nx
import copy
import networker.geo_math as gm
import numpy as np
from rtree import Rtree

"""
Module for GeoGraph extension to networkx Graph class
"""

class GeoObject(object):

    """ 
    Base class for Geo objects 
    
    Attributes:
        srs:  spatial reference system for the coords
            In proj4 string format (http://trac.osgeo.org/proj/)
        coords:  The coordinates of this object 
            This may be a vector or a collection of vectors depending
            on the object type
    
    """

    def __init__(self, srs, coords):
        self.srs = srs
        self.coords = copy.deepcopy(coords)

    def is_geographic(self):
        """ Returns whether the coords are geographic based on proj4 """
        
        sr = osr.SpatialReference()
        sr.ImportFromProj4(self.srs)
        return bool(sr.IsGeographic())


class GeoGraph(GeoObject, nx.Graph):

    """ 
    class representing networkx Graph with Geo components 
    
    Attributes:
        inherited from GeoObject and nx.Graph
        coords:  coords collection associated with node ids by ix
    
    See Also
    --------

    GeoObject
    networkx.Graph
    
    """

    def __init__(self, srs=gm.PROJ4_FLAT_EARTH, coords={}, data=None, **attr):
        """ initialize via both parent classes """
        
        GeoObject.__init__(self, srs, coords)
        nx.Graph.__init__(self, data, **attr)

        # handle case where coords have keys not referenced by edges
        coord_keys = range(len(coords))
        if isinstance(coords, dict):
            coord_keys = coords.keys()
            
        new_nodes = set(coord_keys) - set(self.node.keys())
        self.add_nodes_from(new_nodes)


    def is_aligned(self):
        """
        Test whether this GeoGraph node dict and coords dict are aligned

        Useful in case you have a GeoGraph that has its coords/nodes modified
        """
        
        coord_keys = range(len(self.coords))
        if isinstance(self.coords, dict):
            coord_keys = self.coords.keys()
 
        assert sorted(self.nodes()) == sorted(coord_keys), \
            "GeoGraph nodes and coords not aligned"

        return True

    def project_onto(self, other, rtree_index=None): 
        """
        project other GeoGraph nodes onto their nearest point on this GeoGraph

        Assumes there is no node label overlap between self and other

        Args:
            other:  GeoGraph with nodes to be projected onto this GeoGraph
            rtree_index:  rtree of edges within self to be used for speeding
                up matching of self edges to other nodes

        Returns:
            GeoGraph:  With nearest edges from self to nodes in other 
                and additional nodes (and coords) for projected points on those
                edges (aka 'fake' nodes).  fake nodes will split these edges?
            OR
            dict:  other_node_id -> ((edge), coords) where 
                edge:  (node1, node2) tuple of edge in self other node is 
                    nearest to
                coords:  point on edge other node projects to
        """
        assert len(set(self.nodes()).intersection(set(other.nodes()))) == 0, \
            "the intersection of self and other graphs should be empty"

        projections = {}
        for node in other.nodes():
            edge, coords = self.find_nearest_edge(other.coords[node], rtree_index=rtree_index)
            projections[node] = (edge, coords)
        
        # create new GeoGraph with others coords
        geo = GeoGraph(other.srs, other.coords)

        # add nodes and edges
        # assumes nodes are numeric
        # assign new node ids by sequence starting at max node val of other 
        max_node = max(other.nodes())
        for i, node in enumerate(other.nodes(), 1):
            edge, coords = projections[node]
            new_node = max_node + i 
            # add 'real' node to 'fake' node edge
            geo.add_edge(node, new_node)
            # split edge with new node
            geo.add_edge(edge[0], new_node)
            geo.add_edge(new_node, edge[1])
            # add coords
            geo.coords[new_node] = coords

        return geo 

    def find_nearest_edge(self, coord, rtree_index=None):
        """
        Find the nearest edge to the coordinate in space

        Args:
            coord:  coordinate to lookup nearest edge to
            rtree_index:  rtree of edges within self

        Returns:
            edge:  a tuple of nodes 
            coords:  coordinates of nearest point on edge

        """
        if rtree_index:
            nearest_segment = rtree_index.nearest(np.ravel((coord, coord)), \
                objects=True).next()
            near_edge, coords = nearest_segment.object
            sq_dist, near_coords = self._sq_dist_to_edge(near_edge, coord)
            return near_edge, near_coords

        else:
            # exhaustively compare segments
            min_dist = np.inf
            near_edge = None
            near_coords = None

            for edge in self.edges():
                sq_dist, coords = self._sq_dist_to_edge(edge, coord)
                if sq_dist < min_dist:
                    near_edge = edge
                    min_dist = sq_dist
                    near_coords = coords

            return near_edge, near_coords
     
    def get_coord_edge_set(self):
        """
        get edges as a set of frozensets of coordinate pairs
        """

        tup_map = {i: tuple(coords) for i, coords in self.coords.items()}
        edge_sets = map(lambda e: frozenset([tup_map[e[0]], tup_map[e[1]]]),\
            self.edges())
        return set(edge_sets)


    def _sq_dist_to_edge(self, edge, coord):
        """
        helper to calculate the square distance from coord to
        the edge and the nearest coord on that edge

        Args:
            edge:  tuple of nodes representing edge
            coord:  coordinates of point to take distance to

        Returns:
            square_distance:  from coord to nearest point on edge
            coords:  coordinates of nearest point on edge
        """

        c0 = self.coords[edge[0]]
        c1 = self.coords[edge[1]]

        space = np.shape(coord)[0]
        assert np.shape(c0)[0] == np.shape(c1)[0] == space, \
            "coordinate space of nodes and coord must match"

        project_fun = {2: gm.project_point_on_segment, 
                       3: gm.project_point_on_arc}
        
        proj_coord = project_fun[space](coord, c0, c1)
           
        return np.sum((proj_coord - np.array(coord)) ** 2), proj_coord


    def get_rtree_index(self):
        """
        Get an rtree index of the edges within this GeoGraph
        
        Returns:
            rtree:  rtree of edges indexed by bounding_box (4 coordinates)
                allowing lookup of the edge (node1, node2) and it's segment
        """
        def edge_generator():
            edges_bounds = map(lambda tup: (tup, gm.make_bounding_box(*map(lambda n: self.coords[n], tup))), self.edges())

            for edge, box in edges_bounds:
                # Object is in form of (u.label, v.label), (u.coord, v.coord)
                yield (hash(edge), box, (edge, map(lambda ep: np.array(self.coords[ep]), edge)))


        # Init rtree and store grid edges
        rtree = Rtree(edge_generator())
        return rtree







