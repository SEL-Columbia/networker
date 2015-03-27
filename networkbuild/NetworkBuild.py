# -*- coding: utf-8 -*-

__author__ = 'Brandon Ogle'

from copy import deepcopy

import ogr
import numpy as np
import networkx as nx
import pandas as pd

from rtree import Rtree
from networkbuild.utils import UnionFind

from networkbuild.geo_math import spherical_distance_scalar, \
                                  coordinate_transform_proj4, \
                                  make_bounding_box, \
                                  project_point_to_segment, \
                                  PROJ4_LATLONG

class NetworkBuild(object):

    def __init__(self, metric, existing=None):
        self.existing = existing
        if self.existing is not None:
            grid = self.setup_grid(existing)
            net = self.setup_input_nodes(metric)
            self.G, self.subgraphs, self.rtree = self.grid_settlement_merge(grid, net)
        else:
            self.G = self.setup_input_nodes(metric)

    @staticmethod
    def setup_grid(path):

        #Open the shapefile and get the projection
        driver = ogr.GetDriverByName('ESRI Shapefile')
        shp = driver.Open(path)
        layer = shp.GetLayer()
        spatial_ref = layer.GetSpatialRef()
        input_proj = spatial_ref.ExportToProj4()

        # Load in the grid
        grid = nx.read_shp(path)
        # Convert coord labels to ints
        grid = nx.convert_node_labels_to_integers(grid, label_attribute='coords')

        # if coords in utm convert to latlong
        if input_proj != PROJ4_LATLONG:
            input_coords = np.row_stack(nx.get_node_attributes(grid, 'coords').values())
            coords = coordinate_transform_proj(input_proj, PROJ4_LATLONG, input_coords)
            nx.set_node_attributes(grid, 'coords', {i : coord for i, coord
                                                    in enumerate(coords)})

        # Append 'grid-' to grid node labels
        # TODO:  add 'grid' boolean attribute rather than this

        # grid = nx.relabel_nodes(grid, {n: 'grid-' + str(n) for n in grid.nodes()})
        # Set set grid to True, mv to 0
        nx.set_node_attributes(grid, 'grid', {n:True for n in grid.nodes()})
        nx.set_node_attributes(grid, 'budget', {n:0 for n in grid.nodes()})
        # Set edge weights
        get_coord = lambda x: grid.node[x]['coords']
        nx.set_edge_attributes(grid, 'weight', {(u, v): \
                       spherical_distance_scalar([map(get_coord, [u,v])])  \
                       for u,v in grid.edges()}) 

        return grid.to_undirected()

    @staticmethod
    def setup_input_nodes(path):

        input_proj = csv_projection(path)
        header_row = 1 if input_proj else 0

        # read in the csv
        metrics = pd.read_csv(path, header=header_row)

        # Find the xy pattern used
        coord_pattern = [['x', 'y'], ['X', 'Y'], ['Lon', 'Lat'], ['Longitude', 'Latitude']]
        coord_cols = [[xp, yp] for [xp, yp] in coord_pattern
                if hasattr(metrics, xp) and hasattr(metrics, yp)][0]

        # Stack the coords
        coords = np.column_stack(map(metrics.get, coord_cols))

        #if coords are not in latlong need to convert
        if input_proj and input_proj != PROJ4_LATLONG:
            coords = coordinate_transform_proj(input_proj, PROJ4_LATLONG, coords)

        # Store the MV array
        if hasattr(metrics, 'Demand > Projected nodal demand per year'):
            mv = metrics['Demand > Projected nodal demand per year'].values
        else:
            mv = metrics['Demand...Projected.nodal.demand.per.year'].values

        # Setup the network
        net = nx.Graph()
        net.add_nodes_from(range(mv.size))

        # Store the metric attrs
        nx.set_node_attributes(net, 'coords', {n: coords[n] for n in range(mv.size)})
        nx.set_node_attributes(net, 'budget', {n: mv[n] for n in range(mv.size)})

        return net

    @staticmethod
    def project_to_grid(rtree, coords):
        """
        Reduces the point set to sqrt(n) via kmeans, 
        then projects those representative points onto the grid
        """
        fake_nodes = []
        # neighborhoods = kmeans(coords, np.sqrt(coords.shape[0]))[0]
        # for coord in neighborhoods:
        for coord in coords:
            # Find nearest bounding box
            nearest_segment = rtree.nearest(np.ravel((coord, coord)), objects=True).next()
            uv, line = nearest_segment.object
            # project the coord onto the edge and append to the list of fake nodes
            fake = project_point_to_segment(coord, *line)
            #Todo: only add unique fake nodes
            fake_nodes.append((uv, fake))

        return fake_nodes

    @staticmethod
    def grid_settlement_merge(grid, net):
        """
        Sets up the Graph, UnionFind (DisjoinSet), and RTree datastructures
        for use in network algorithms

        Args:
            grid:  graph representing existing grid (assumes node ids don't conflict
                with net (demand) nodes)
            net:  graph of nodes representing demand

        Returns:
            graph:  graph with demand nodes and their nearest nodes to the 
                existing grid (i.e. 'fake' nodes)
            subgraphs:  UnionFind datastructure populated with fake nodes and
                associated with the appropriate connected component
            rtree:  spatial index populated with the edges from the 
                existing grid
        """

        # Coordinates of grid vertices, lookup strips the 'grid-' prior to indexing
        g_coords = nx.get_node_attributes(grid, 'coords')

        def edge_generator():
            edges_bounds = map(lambda tup: (tup, make_bounding_box(*map(lambda n: g_coords[n], tup))), grid.edges())

            for edge, box in edges_bounds:
                # Object is in form of (u.label, v.label), (u.coord, v.coord)
                yield (hash(edge), box, (edge, map(lambda ep: np.array(g_coords[ep]), edge)))


        # Init rtree and store grid edges
        rtree = Rtree(edge_generator())

        # Project Fake Nodes
        coords = np.row_stack(nx.get_node_attributes(net, 'coords').values())
        fake_nodes = NetworkBuild.project_to_grid(rtree, coords)

        # Get the grid components to init mv grid centers
        subgrids = nx.connected_components(grid)

        # Init the DisjointSet
        subgraphs = UnionFind()

        # Build the subgrid components
        for sub in subgrids:
            # Start subcomponent with first node
            subgraphs.add_component(sub[0], budget=grid.node[sub[0]]['budget'])
            # Merge remaining nodes with component
            for node in sub[1:]:
                subgraphs.add_component(node, budget=grid.node[node]['budget'])
                # The existing grid nodes have are on the grid 
                # (so distance is 0)
                subgraphs.union(sub[0], node, 0)

        # setup graph to be populated with fake nodes
        big = deepcopy(net)
        
        # Id of the last real node
        last_real = max(net.nodes())

        # Fake_nodes
        for idx, ((u, v), fake) in enumerate(fake_nodes, 1):

            # Make sure something wonky isn't going on
            assert(subgraphs[u] == subgraphs[v])

            # Add the fake node to the big net
            fake_id = last_real + idx
            big.add_node(fake_id, coords=fake, budget=np.inf)

            # Merge the fake node with the grid subgraph
            subgraphs.add_component(fake_id, budget=np.inf)
            subgraphs.union(fake_id, u, 0)

        return big, subgraphs, rtree


    def build(self, func, min_node_component=2):
        if self.existing:
            opt_net = func(self.G, self.subgraphs, self.rtree)
        else:
            opt_net = func(self.G)


        return nx.union_all(filter(lambda sub: len(sub.node) >= min_node_component,
            nx.connected_component_subgraphs(opt_net)))
