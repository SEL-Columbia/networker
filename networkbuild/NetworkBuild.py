# -*- coding: utf-8 -*-
__author__ = 'Brandon Ogle'

import numpy as np
import networkx as nx
import pandas as pd

from rtree import Rtree
from shapely.geometry import LineString, Point
from networkbuild.utils import UnionFind, make_bounding_box

class NetworkBuild(object):

    def __init__(self, metric, existing=None):
        self.existing = existing
        if self.existing is not None:
            grid = self.setup_grid(existing)
            net  = self.setup_input_nodes(metric)
            self.G, self.subgraphs, self.rtree = self.grid_settlement_merge(grid, net)
        else:
            self.G = self.setup_input_nodes(metric)
    
    @staticmethod
    def setup_grid(path):
    
        # Load in the grid
        grid = nx.read_shp(path)
        # Convert coord labels to ints
        grid = nx.convert_node_labels_to_integers(grid, label_attribute='coords')
        # Append 'grid-' to grid node labels
        grid = nx.relabel_nodes(grid, {n: 'grid-' + str(n) for n in grid.nodes()})
        # Set mv to 0
        nx.set_node_attributes(grid, 'mv', {n:0 for n in grid.nodes()})
        
        return grid.to_undirected()

    @staticmethod
    def setup_input_nodes(path):
        
        # Load in the metrics file
        metrics = pd.read_csv(path)
        # Stack the coords
        coords = np.column_stack((metrics.x, metrics.y))
        # Store the MV array
        mv = metrics['Demand...Projected.nodal.demand.per.year'].values
        
        # Setup the network
        net = nx.Graph()
        net.add_nodes_from(range(mv.size))
        
        # Store the metric attrs
        nx.set_node_attributes(net, 'coords', {n: coords[n] for n in range(mv.size)})
        nx.set_node_attributes(net, 'mv', {n: mv[n] for n in range(mv.size)})
        
        return net

    @staticmethod
    def grid_settlement_merge(grid, net):
        
        # Coordinates of grid vertices, lookup strips the 'grid-' prior to indexing
        g_coords_array = np.row_stack(nx.get_node_attributes(grid, 'coords').values())
        g_coords = lambda x: g_coords_array[int(x.strip('grid-'))]
        edges_bounds = map(lambda tup: (tup, make_bounding_box(*map(lambda n: g_coords(n), tup))), grid.edges())
        
        # init rtree and store grid edges
        rtree = Rtree()
        for edge, box in edges_bounds:
            # Object is in form of (u.label, v.label), (u.coord, v.coord)
            rtree.insert(hash(edge), box, obj=(edge, map(g_coords, edge)))
        
        # project fake nodes
        fake_nodes = []
        coords = np.row_stack(nx.get_node_attributes(net, 'coords').values())
        
        for coord in coords:
            # Find nearest bounding box
            nearest_segment = rtree.nearest(np.ravel((coord, coord)), objects=True).next()
            uv, line = nearest_segment.object
            # project the coord onto the edge and append to the list of fake nodes
            edge = LineString(line)
            fake = edge.interpolate(edge.project(Point(coord)))
            fake_nodes.append((uv, np.asarray(fake.coords.xy).reshape(2)))
        
        # Get the grid components to init mv grid centers
        subgrids = nx.connected_components(grid)
        
        # Union the subgraphs, and init the DisjointSet
        big = nx.union(grid, net)
        subgraphs = UnionFind(big)
        
        # Build the subgrid components
        for sub in subgrids:
            # Start subcomponent with first node
            subgraphs[sub[0]]
            # Merge remaining nodes with component
            for node in sub[1:]:
                subgraphs[node]
                # The existing grid nodes have dummy mv
                subgraphs.union(sub[0], node, 0)

        # Id of the last real node
        last_real = max(net.nodes()) 
        
        # Fake_nodes
        for idx, ((u, v), fake) in enumerate(fake_nodes, 1):

            # Make sure something wonky isn't going on
            assert(subgraphs[u] == subgraphs[v])
                
            # Add the fake node to the big net
            fake_id = last_real + idx
            big.add_node(fake_id, coords=fake, mv=np.inf)
            
            # Merge the fake node with the grid subgraph
            subgraphs.union(fake_id, u, 0)
            
        # Validate output
        
        # assert the number of disjoint sets == subcomponents with more than one node
        assert(len(subgraphs.connected_components()) == len(filter(lambda x: len(x) > 1, nx.connected_components(big))))    
        
        # Now that the needed grid data has been captured
        # Pull the existing grid out of the network
        big.remove_edges_from(grid.edges())
        big.remove_nodes_from(grid.nodes())
        
        return big, subgraphs, rtree

    
    def build(self, func, min_node_component=2):
        if self.existing:
            opt_net = func(self.G, self.subgraphs, self.rtree)
        else:
            opt_net = func(self.G)

        
        return nx.union_all(filter(lambda sub: len(sub.node) >= min_node_component, 
            nx.connected_component_subgraphs(opt_net)))

