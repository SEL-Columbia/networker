# -*- coding: utf8 -*-
__author__ = 'Brandon Ogle'

import heapq
import osr
import numpy as np
import networkx as nx

from collections import defaultdict
from numba import jit

class Dict(dict):
    """dictionary allowing weakref"""
    pass

class UnionFind:

    def __init__(self):
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
        self.weights = {}
        self.budget = {}
        self.parents = {}
        self.children = Dict() #This was previously used such that modifying it changed all refs pointing here
        self.queues = {}
    
    def __getitem__(self, object):
        """Find and return the name of the set containing the object."""
        # check for previously unknown object
        if object not in self.parents:
            self.add_component(object)
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

    def add_component(self, object, budget=0, weight=1):
        """ 
        Add standalone component to the UnionFind structure by 
        initializing it's dict values.

        Note:  Need to use this method to specify budget

        Args:
            object:  identifier for the component
            budget:  The budget allocated to this node for connection costs
            weight:  The weight associated with this node for agglomeration 
        """
        # only add if not already here
        if object not in self.parents:
            self.parents[object] = object
            self.weights[object] = weight 
            self.budget[object] = budget
            self.children[object] = [object]
            self.queues[object] = PriorityQueue()
 

    def __iter__(self):
        """Iterate through all items ever found or unioned by this structure."""
        return iter(self.parents)
    

    def push(self, queue, item, priority):
        """Pushes an item into component queue"""
        u, v = item
        queue.push(item, priority)


    def union(self, g1, g2, d):
        """
        Find the sets containing the objects and merge them all.
        During merge; children, budget and queues are all aggregated
        into the dominant or 'heaviest' set.

        Args:
            g1 (obj): Key of a member in the disjoint sets
            g2 (obj): Key of a member in the disjoint sets
             d (Num): distance between g1 and g2
        """
        graphs = [g1, g2]; 
        point_budget = np.array(map(self.budget.get, graphs))
        fake = point_budget == np.inf

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
            # if the fake node is also a parent component
            if grid == smallest:
                # only set the parent budget, leave fake node budget alone so that
                # it continues to have 'infinite' budget
                self.budget[heaviest] -= d
                return

        self.budget[heaviest] += self.budget[smallest] - d
        self.budget[smallest] = self.budget[heaviest]


    def connected_components(self, component_subset=None):
        """Return the roots for all disjoint sets
        
        Args:
            component_subset:  a set of nodes for which all connected_component
              children must be within (if None, then return all)

        Returns:
            connected_components
        """
        if component_subset:
            return set([self[r] for r in self if \
                any(c in component_subset for c in self.children[self[r]])])
        else:
            return set([self[r] for r in self])


    def component_set(self, component):
        """Return the component set of the objects

        Args:
            n (obj): member node in one of the disjoint sets

        Returns:
            List of nodes in the same set as n
        """

        return [c for c in self.children[self[component]]]


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

# Utility functions
def csv_projection(path):
    """
    Get projection from csv file of geo data

    Args: 
        path (str): path to the csv

    Returns:
        Proj4 projection if included in header else None
    """

    with open(path) as raw_text:
        header = raw_text.readline()

    if 'PROJ.4' in header:
        return header

# plot maps
def draw_np_graph(g, node_color='r', edge_color='b', node_label_field='ix', edge_label_field=None, node_size=300):

    from mpl_toolkits.basemap import Basemap
    m = Basemap(
            projection='merc',
            ellps = 'WGS84',
            llcrnrlon=0,
            llcrnrlat=0,
            urcrnrlon=1,
            urcrnrlat=1,
            lat_ts=0,
            resolution='i',
            suppress_ticks=True)

    node_pos = {nd: m(g.node[nd]['coords'][0], g.node[nd]['coords'][1]) for nd in g.nodes()}
        
    # handle main graph rendering 
    if node_label_field:
        if node_label_field != 'ix':
            label_vals = nx.get_node_attributes(g, node_label_field)
            nx.draw_networkx(g, pos=node_pos, labels=label_vals, node_color=node_color, edge_color=edge_color, node_size=node_size)

        else: # label via ids
            nx.draw_networkx(g, pos=node_pos, with_labels=True, node_color=node_color, edge_color=edge_color, node_size=node_size)

    else:
        nx.draw_networkx(g, pos=node_pos, with_labels=False, node_color=node_color, edge_color=edge_color, node_size=node_size)

    # handle edge labels if needed
    if edge_label_field:
        edge_labels = nx.get_edge_attributes(g, edge_label_field)
        nx.draw_networkx_edge_labels(g, pos=node_pos, edge_labels=edge_labels)


# to js
def network_to_json(g):

    from mpl_toolkits.basemap import Basemap
    from networkx.readwrite import json_graph
    m = Basemap(
            projection='merc',
            ellps = 'WGS84',
            llcrnrlon=0,
            llcrnrlat=0,
            urcrnrlon=1,
            urcrnrlat=1,
            lat_ts=0,
            resolution='i',
            suppress_ticks=True)

    g2 = g.copy()
    for nd in g.nodes():
        g2.node[nd]['coords'] = m(g.node[nd]['coords'][0], g.node[nd]['coords'][1])

    js_g = json_graph.node_link_data(g2)
    return js_g
 
