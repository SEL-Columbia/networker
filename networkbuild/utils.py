# -*- coding: utf8 -*-
__author__ = 'Brandon Ogle'

import heapq

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
        
        self.parents[r] = heaviest
        self.weights[heaviest] += self.weights[r]
        
        self.mv[heaviest] = (self.mv[g1] + self.mv[g2]) - d
        self.mv[r] = self.mv[heaviest]
        
        self.children[heaviest] += self.children[r]
        self.children[r] = self.children[heaviest]
    
        self.queues[heaviest].merge(self.queues[r])
        self.queues[r] = self.queues[heaviest] 
    
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
