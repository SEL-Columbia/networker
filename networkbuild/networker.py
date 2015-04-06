# -*- coding: utf-8 -*-

import networkx as nx
import pandas as pd

from networkbuild.utils import UnionFind, csv_projection
from networkbuild.classes import GeoGraph
import networkbuild.network_io as nio
import networkbuild.geo_math as gm
import networkbuild.algorithms as algo

class Networker(object):

    """
    class for running minimum spanning forest algorithms on a set of
    spatially referenced nodes

    Attributes:
        config:  dict (potentially nested) of configuration params
            TODO:  detail the params and define defaults 

            existing_networks:  
                filename:  shapefile
                budget_value:  { 0 }
            demand_node_metrics:  
                filename:  
                x_column: column for x or longitude values { 'X' }
                y_column: column for y or latitude values { 'Y' } 
                budget_column:  column for nodal budget { 'metric' }
            network_algorithm: mod_boruvka, mod_kruskal
            network_parameters:
                minimum_node_count:  min number of nodes in sub-network { 2 } 
            output_directory: { 'output' }
    """

    ALGOS = {'mod_boruvka': algo.mod_boruvka}

    def __init__(self, config):
        self.config = config

    def run(self):
        """
        run a minimum spanning forest algorithm on inputs and write output 
        based on configuration
        """

        geo_graph = subgraphs = rtree = None

        demand_nodes = load_node_metrics(**config['demand_node_metrics'])
        if config.has_key('existing_networks'): 
            existing = load_existing_networks(**config['existing_networks'])
            # rename existing nodes so that they don't intersect with metrics 
            existing = nx.relabel_nodes(existing, \
                {n: 'grid-' + str(n) for n in existing.nodes()})
            geo_graph, subgraphs, rtree = merge_network_and_nodes(existing, demand_nodes)
        else:
            geo_graph = demand_nodes

        # now run the selected algorithm
        network_algo = ALGOS[config['network_algorithm']]
        result_geo_graph = network_algo(geo_graph, subgraphs=subgraphs, rtree=tree)

        # filter out any unused fake nodes

        # now filter out subnetworks via minimum node count
        min_node_count = config['network_parameters']['minimum_node_count']
        filtered_graph = nx.union_all(filter(lambda sub: len(sub.node) >= min_node_count,
            nx.connected_component_subgraphs(result_geo_graph)))
        # map coords back to geograph
        coords = {i: result_geo_graph.coords[i] for i in filtered_graph}
        geo = GeoGraph(result_geo_graph.srs, coords=coords, data=filtered_graph)

        # now save it
        nio.write_shp(geo, config['output_directory'])
            
    @staticmethod
    def merge_network_and_nodes(network, demand_nodes):
        """
        merge the network and nodes GeoGraphs to set up the Graph, UnionFind 
        (DisjoinSet), and RTree datastructures for use in network algorithms

        Args:
            network:  graph representing existing network (assumes node ids don't conflict
                with net (demand) nodes)
            demand_nodes:  graph of nodes representing demand

        Returns:
            graph:  graph with demand nodes and their nearest nodes to the 
                existing network (i.e. 'fake' nodes)
            subgraphs:  UnionFind datastructure populated with fake nodes and
                associated with the appropriate connected component from the
                network
            rtree:  spatial index populated with the edges from the 
                existing network 
        """
        
        # project demand nodes onto network
        rtree = network.get_rtree_index()
        grid_with_fakes = network.project_onto(demand_nodes, rtree_index=rtree)
        
        # get only the fake nodes and the associated network edges
        demand_node_set = set(demand_nodes.nodes())
        net_plus_demand = set(network.nodes()) + demand_node_set
        fakes = set(grid_with_fakes.nodes()) - net_plus_demand_
        # fake node should only have 2 neighbors from the existing network
        # that is the nearest edge
        get_fake_edge = lambda node: tuple(set(grid_with_fakes.neighbors(x)) - \
            demand_node_set)
        edge_fakes = [(get_fake_edge(fake), fake) for fake in fakes]

        # Get the network components to init budget centers
        subnets = nx.connected_components(network)

        # Init the DisjointSet
        subgraphs = UnionFind()

        # Build the subnet components
        for sub in subnets:
            # Start subcomponent with first node
            subgraphs.add_component(sub[0], budget=network.node[sub[0]]['budget'])
            # Merge remaining nodes with component
            for node in sub[1:]:
                subgraphs.add_component(node, budget=network.node[node]['budget'])
                # The existing grid nodes are on the grid (so distance is 0)
                subgraphs.union(sub[0], node, 0)

        # setup merged graph to be populated with fake nodes
        merged = GeoGraph(demand_nodes.srs, demand_nodes.coords)
        # merge fakes in
        for ((u, v), fake) in edge_fakes:

            # Make sure something wonky isn't going on
            assert(subgraphs[u] == subgraphs[v])

            # Add the fake node to the big net
            merged.add_node(fake, budget=np.inf)
            merged.coords[fake] = grid_with_fakes.coords[fake]

            # Merge the fake node with the grid subgraph
            subgraphs.add_component(fake, budget=np.inf)
            subgraphs.union(fake, u, 0)

        return merged, subgraphs, rtree


    @staticmethod
    def load_existing_networks(filename="existing_networks.shp", budget_value=0):
        """
        load existing_networks shp into GeoGraph nodes, edges and assign budget 

        Args:
            filename:  existing_networks shapefile
            budget_value:  default budget value for nodes in existing network
        
        Returns:
            GeoGraph of existing networks with budget attribute
        """
 
        geo_net = nio.load_shp(filename)
        
        nx.set_node_attributes(geo_net, 'budget', {n:budget_value for n in geo_net.nodes()})

        # determine edge weight/distance function by whether they're geocentric
        # (assuming coordinates are not in 3-space here)
         
        distance = gm.spherical_distance if geo_net.is_geographic else \
            gm.euclidean_distance 

        nx.set_edge_attributes(geo_net, 'weight', {(u, v): \
                       distance(map(geo_net.coords.get, [u,v]))  \
                       for u,v in geo_net.edges()}) 

        return geo_net


    @staticmethod
    def load_node_metrics(filename="metrics.csv", x_column="X", y_column="Y", \
        budget_column="metric"):
        """
        load node_metrics csv into GeoGraph (nodes only)

        Args:
            filename:  nodal metrics csv file
            x_column, y_column:  col names to take x, y from
            budget_column:  col to take budget from
        
        Returns:
            GeoGraph of nodes only with budget attribute
        """
    
        input_proj = csv_projection(metric_config['file'])
        header_row = 1 if input_proj else 0

        # read in the csv
        metrics = pd.read_csv(path, header=header_row)

        coord_cols = [x_column, y_column]
        assert all([hasattr(metrics, col) for col in coord_cols]), \
            "metrics file does not contain coordinate columns {}, {}".\
            format(x_column, y_column)

        assert hasattr(metrics, budget_column), \
            "metrics file does not contain budget columns {}".\
            format(budget_column)

        # Stack the coords and get the budget
        coords = np.column_stack(map(metrics.get, coord_cols))
        budget = metrics[budget_column]

        to_dict = lambda column: {i, value for i, value in enumerate(column)}
        coords_dict = to_dict(coords)
        budget_dict = to_dict{budget}

        geo_nodes = GeoGraph(input_proj, coords_dict)
        
        nx.set_node_attributes(geo_nnodes, 'budget', budget_dict)
        return geo_nodes
       

    @staticmethod
    def _validate(config):
        """
        validate configuration
        """
