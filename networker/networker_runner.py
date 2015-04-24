# -*- coding: utf-8 -*-

import networkx as nx
import numpy as np
import pandas as pd

from networker.classes.unionfind import UnionFind
from networker.utils import csv_projection
from networker.classes.geograph import GeoGraph
import networker.io as nio
import networker.geomath as gm
import networker.algorithms as algo
import os
import json
import jsonschema


class NetworkerRunner(object):

    """
    class for running minimum spanning forest algorithms on a set of
    spatially referenced nodes

    Attributes:
        config:  nested dict of configuration params

            existing_networks:
                filename:  shapefile
                budget_value:  { 0 }
            demand_nodes:
                filename:
                x_column: column for x or longitude values { 'X' }
                y_column: column for y or latitude values { 'Y' }
                budget_column:  column for nodal budget { 'metric' }
            network_algorithm: mod_boruvka, mod_kruskal
            network_parameters:
                minimum_node_count:  min number of nodes in sub-network { 2 }
    """

    ALGOS = {'mod_boruvka': algo.mod_boruvka}
    SCHEMA_FILE = "networker_config_schema.json"

    def __init__(self, config, output_directory="."):
        self.config = config
        self.output_directory = output_directory

    def run(self):
        """
        run a minimum spanning forest algorithm on inputs and write output
        based on configuration
        """

        demand_nodes = load_node_metrics(**self.config['demand_nodes'])
        existing_networks = None
        if 'existing_networks' in self.config:
            existing_networks = load_existing_networks(
                prefix="grid-",
                **self.config['existing_networks'])

        min_node_count = self.config['network_parameters']\
                                    ['minimum_node_count']
        network_algorithm = self.config['network_algorithm']

        msf, _ = build_network(demand_nodes, 
                                existing=existing_networks,
                                min_node_count=min_node_count,
                                network_algorithm=network_algorithm)

        # now save it
        if not os.path.exists(self.output_directory):
            os.makedirs(self.output_directory)

        nio.write_shp(msf, self.output_directory)


    def validate(self):
        """
        validate configuration
        throws jsonschema Validate exception if invalid
        """

        # load schema and validate it via jsonschema
        schema_path = os.path.join(os.path.dirname(
            os.path.abspath(__file__)), NetworkerRunner.SCHEMA_FILE)
        schema = json.load(open(schema_path))
        jsonschema.validate(self.config, schema)


def build_network(demand_nodes, 
                    existing=None, 
                    min_node_count=2,
                    network_algorithm='mod_boruvka',
                    one_based=False 
                    ):
    """
    project demand nodes onto optional existing supply network and
    return the 'optimized' network

    Args:
        demand_nodes:  GeoGraph of demand nodes
        existing:  GeoGraph of existing grid (assumes node ids
            don't conflict with demand_nodes
        min_node_count:  minimum number of nodes allowed in a subgraph
            of the result
        network_algorithm:  Algorithm from ALGOS to run
        one_based:  Whether result GeoGraph's nodes should be one_based
            (if not, they are 0 based)

    Returns:
        msf: GeoGraph of minimum spanning forest proposed by the chosen
            network algorithm
        existing: The existing grid GeoGraph (None if it doesn't exist) 
        
    """
    geo_graph = subgraphs = rtree = None

    if existing:
        geo_graph, subgraphs, rtree = \
            merge_network_and_nodes(existing, demand_nodes)
    else:
        geo_graph = demand_nodes

    # now run the selected algorithm
    network_algo = NetworkerRunner.ALGOS[network_algorithm]
    result_geo_graph = network_algo(geo_graph, subgraphs=subgraphs,
                                    rtree=rtree)

    # TODO: Remove unreferenced fake nodes?

    # now filter out subnetworks via minimum node count
    # TODO:  update union_all to support GeoGraph?
    filtered_graph = nx.union_all(filter(
        lambda sub: len(sub.node) >= min_node_count,
        nx.connected_component_subgraphs(result_geo_graph)))

    # map coords back to geograph
    # NOTE:  explicit relabel to int as somewhere in filtering above, some
    # node ids are set to numpy types which screws up comparisons to tuples
    # in write op
    # NOTE:  relabeling nodes in-place here drops node attributes for some
    #   reason so create a copy for now
    def id_label(i): 
        id = int(i+1) if one_based else int(i)
        return id

    msf = None
    if filtered_graph:
        coords = {id_label(i): result_geo_graph.coords[i]
                    for i in filtered_graph}
        relabeled = nx.relabel_nodes(filtered_graph, {i: id_label(i)
            for i in filtered_graph}, copy=True)
        msf = GeoGraph(result_geo_graph.srs, coords=coords, data=relabeled)

    return msf, existing


def merge_network_and_nodes(network, demand_nodes):
    """
    merge the network and nodes GeoGraphs to set up the Graph, UnionFind
    (DisjoinSet), and RTree datastructures for use in network algorithms

    Args:
        network:  graph representing existing network
            (assumes node ids don't conflict with net (demand) nodes)
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
    net_plus_demand = set(network.nodes()).union(demand_node_set)
    fakes = set(grid_with_fakes.nodes()) - net_plus_demand
    # fake node should only have 2 neighbors from the existing network
    # that is the nearest edge
    def get_fake_edge(node): return tuple(set(
        grid_with_fakes.neighbors(node)) - demand_node_set)

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
    merged = GeoGraph(demand_nodes.srs, demand_nodes.coords, data=demand_nodes)
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


def load_existing_networks(filename="existing_networks.shp", budget_value=0,
                           prefix=None):
    """
    load existing_networks shp into GeoGraph nodes, edges and assign budget

    Args:
        filename:  existing_networks shapefile
        budget_value:  default budget value for nodes in existing network
        prefix: if not None, relabel node ids with the prefix

    Returns:
        GeoGraph of existing networks with budget attribute
    """

    geo_net = nio.load_shp(filename)

    nx.set_node_attributes(geo_net, 'budget', {n: budget_value
        for n in geo_net.nodes()})

    # determine edge weight/distance function by whether they're geocentric
    # (assuming coordinates are not in 3-space here)

    distance = gm.spherical_distance if geo_net.is_geographic else \
        gm.euclidean_distance

    nx.set_edge_attributes(geo_net, 'weight', {(u, v):
                   distance(map(geo_net.coords.get, [u, v]))
                   for u, v in geo_net.edges()})

    if prefix:
        nx.relabel_nodes(geo_net,
            {n: prefix + str(n) for n in geo_net.nodes()}, copy=False)
        geo_net.coords = {prefix + str(n): c for n, c in
            geo_net.coords.items()}

    return geo_net


def load_node_metrics(filename="metrics.csv", x_column="X", y_column="Y",
                        budget_column="metric", budget_value=1000):
    """
    load node_metrics csv into GeoGraph (nodes only)

    Args:
        filename:  nodal metrics csv file
        x_column, y_column:  col names to take x, y from
        budget_column:  col to take budget from
        budget_value: default value for nodal budget

    Returns:
        GeoGraph of nodes only with budget attribute
    """

    input_proj = csv_projection(filename)
    header_row = 1 if input_proj else 0

    # read in the csv
    # NOTE:  Convert x,y via float cast to preserve precision of input
    metrics = pd.read_csv(filename, header=header_row,
        converters={x_column: float, y_column: float})

    coord_cols = [x_column, y_column]
    assert all([hasattr(metrics, col) for col in coord_cols]), \
        "metrics file does not contain coordinate columns {}, {}".\
        format(x_column, y_column)

    # default budget
    budget = len(metrics) * [budget_value]
    # try to get nodal budgets
    if hasattr(metrics, budget_column):
        budget = metrics[budget_column].fillna(budget_value)

    # Stack the coords and get the budget
    coords = np.column_stack(map(metrics.get, coord_cols))

    # set default projection
    if not input_proj:
        if gm.is_in_lon_lat(coords):
            input_proj = gm.PROJ4_LATLONG
        else:
            input_proj = gm.PROJ4_FLAT_EARTH

    coords_dict = dict(enumerate(coords))
    budget_dict = dict(enumerate(budget))

    geo_nodes = GeoGraph(input_proj, coords_dict)

    nx.set_node_attributes(geo_nodes, 'budget', budget_dict)
    return geo_nodes
