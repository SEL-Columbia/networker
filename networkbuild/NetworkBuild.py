# -*- coding: utf-8 -*-
__author__ = 'Brandon Ogle'

import ogr
import numpy as np
import networkx as nx
import pandas as pd

from rtree import Rtree
from networkbuild.utils import UnionFind, make_bounding_box, project_point,\
                               csv_projection, string_to_proj4, utm_to_wgs84

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
        spatialRef = layer.GetSpatialRef()
        proj4 = string_to_proj4(spatialRef.ExportToProj4())

        # Load in the grid
        grid = nx.read_shp(path)
        # Convert coord labels to ints
        grid = nx.convert_node_labels_to_integers(grid, label_attribute='coords')

        # if coords in utm convert to latlong
        if proj4['proj'] == 'utm':
            utmcoords = np.row_stack(nx.get_node_attributes(grid, 'coords').values())
            coords = utm_to_wgs84(utmcoords, int(proj4['zone']))
            nx.set_node_attributes(grid, 'coords', {i : coord for i, coord
                                                    in enumerate(coords)})

        # Append 'grid-' to grid node labels
        grid = nx.relabel_nodes(grid, {n: 'grid-' + str(n) for n in grid.nodes()})
        # Set mv to 0
        nx.set_node_attributes(grid, 'mv', {n:0 for n in grid.nodes()})

        return grid.to_undirected()

    @staticmethod
    def setup_input_nodes(path):

        proj4 = csv_projection(path)
        header_row = 1 if proj4 else 0

        # read in the csv
        metrics = pd.read_csv(path, header=header_row)

        # Find the xy pattern used
        coord_pattern = [['x', 'y'], ['X', 'Y'], ['Lon', 'Lat'], ['Longitude', 'Latitude']]
        coord_cols = [[xp, yp] for [xp, yp] in coord_pattern
                if hasattr(metrics, xp) and hasattr(metrics, yp)][0]

        # Stack the coords
        coords = np.column_stack(map(metrics.get, coord_cols))

        #if coords are in utm, need to convert to longlat
        if proj4['proj'] == 'utm':
            coords = utm_to_wgs84(coords, int(proj4['zone']))

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
        nx.set_node_attributes(net, 'mv', {n: mv[n] for n in range(mv.size)})

        return net

    @staticmethod
    def project_to_grid(rtree, coords):
        """
        Naively projects all nodes in the graph
        onto the existing grid.
        """
        fake_nodes = []

        for coord in coords:
            # Find nearest bounding box
            nearest_segment = rtree.nearest(np.ravel((coord, coord)), objects=True).next()
            uv, line = nearest_segment.object
            # project the coord onto the edge and append to the list of fake nodes
            fake = project_point(line, coord)
            #Todo: only add unique fake nodes
            fake_nodes.append((uv, fake))

        return fake_nodes

    @staticmethod
    def grid_settlement_merge(grid, net):
        # Coordinates of grid vertices, lookup strips the 'grid-' prior to indexing
        g_coords = nx.get_node_attributes(grid, 'coords')

        def edge_generator():
            edges_bounds = map(lambda tup: (tup, make_bounding_box(*map(lambda n: g_coords[n], tup))), grid.edges())

            for edge, box in edges_bounds:
                # Object is in form of (u.label, v.label), (u.coord, u.coord)
                yield (hash(edge), box, (edge, map(lambda ep: np.array(g_coords[ep]), edge)))


        # Init rtree and store grid edges
        rtree = Rtree(edge_generator())

        # Project Fake Nodes
        coords = np.row_stack(nx.get_node_attributes(net, 'coords').values())
        fake_nodes = NetworkBuild.project_to_grid(rtree, coords)

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
