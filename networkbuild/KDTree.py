__author__ = 'Brandon Ogle'

import numpy as np

def distance(u, v): 
    return np.sum((u - v)**2)

class KDTree(object):    
    def __init__(self, data, depth=0, make_idx=True):
        self.n, self.k = data.shape
        
        if make_idx:
            # index the data
            data = np.column_stack((data, np.arange(self.n)))
        else:
            # subtract the indexed dimension in later calls
            self.k -= 1
        
        self.build(data, depth)
        
    def build(self, data, depth):
        
        if data.size > 0:
            # get the axis to pivot on
            self.axis = depth % self.k
            # sort the data
            s_data = data[np.argsort(data[:, self.axis])]
            # find the pivot point
            point = s_data[len(s_data) // 2]
            
            # point coord
            self.point = point[:-1]
            # point index
            self.idx = int(point[-1])
            
            # all nodes below this node
            self.children = s_data[np.all(s_data[:, :-1] != self.point, axis=1)]
            # branches
            self.left  = KDTree(s_data[: len(s_data) // 2], depth+1, False)
            self.right = KDTree(s_data[len(s_data) // 2 + 1: ], depth+1, False)
        else:
            # empty node
            self.axis=0
            self.idx = self.point = self.idx = self.left = self.right = None
            self.children = np.array([])
            
    def query(self, point, best=None):

        if self.point is None:
            return best
        
        # Dead end backtrack up the tree
        if best is None:
            best = (self.idx, self.point)
        
        # check if current node is closer than best
        if distance(self.point, point) < distance(best[1], point):
            best = (self.idx, self.point)

        # continue traversing the tree
        best = self.near_tree(point).query(point, best)
          
        # traverse the away branch if the orthogonal distance is less than best
        if self.orthongonal_dist(point) < distance(best[1], point):
            best = self.away_tree(point).query(point, best)    
        return best 
    
    def query_subset(self, point, subset, best=None):
        
        # if point in subset, update best
        if self.idx in subset:
            # if closer than current best, or best is none update
            if best is None or distance(self.point, point) < distance(best[1], point):
                best = (self.idx, self.point)
        
        # Dead end backtrack up the tree
        if self.point is None:
            return best
        
        near = self.near_tree(point)
        far = self.away_tree(point)
        
        # what nodes are in the near branch
        if near.children.size > 1:
            near_set = {near.idx}.union(near.children[:, -1])   
        else: near_set = {near.idx}
            
        # check the near branch, if its nodes intersect with the queried subset
        # otherwise move to the away branch
        if set(subset).intersection(near_set):
            best = near.query_subset(point, subset, best)
        else:
            best = far.query_subset(point, subset, best)
        
        # validate best, by ensuring closer point doesn't exist just beyond partition
        if best is not None:
            if self.orthongonal_dist(point) < distance(best[1], point):
                best = far.query_subset(point, subset, best)    
        
        return best 
             
    def orthongonal_dist(self, point):
        orth_point = np.copy(point)
        orth_point[self.axis] = self.point[self.axis]
        return distance(point, self.point)
        
    def near_tree(self, point):
        if point[self.axis] < self.point[self.axis]:
            return self.left
        return self.right
        
    def away_tree(self, point):
        if self.near_tree(point) == self.left:
            return self.right
        return self.left
