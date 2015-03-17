# -*- coding: utf-8 -*-
__author__ = 'Brandon Ogle'

import json
import os, sys
import networkx as nx
import numpy as np

from np.lib import dataset_store, metric, network, variable_store as VS
from networkbuild import NetworkBuild, modBoruvka
np.warnings.filterwarnings('ignore')

class NP_Boruvka_Interface(object):

    def __init__(self, metrics, metric_conf, output, grid=None):

        self.output = output
        # make output dir if not exists
        if not os.path.exists(self.output):
            os.makedirs(self.output)

        targetPath = os.path.join(self.output, "database.db")
        self.dstore = dataset_store.create(targetPath, metrics)

        # setup models
        self.metricModel = metric.getModel('mvMax5')
        self.metricConfiguration = json.load(open(metric_conf))
        self.existing_grid = grid

        self.metricValueByOptionBySection = self.dstore.applyMetric(self.metricModel,
                                                                    self.metricConfiguration)

        self.build_network()
        self.update_metrics()
        self.save_output()

    def update_metrics(self):
        self.metricValueByOptionBySection = self.dstore.updateMetric(self.metricModel,
                                                                     self.metricValueByOptionBySection)

    def build_network(self):
        # load in the grid
       
        grid = nx.Graph() # empty existing grid by default
        G = self._metric_graph()
        DS = R = None
        if(self.existing_grid):
            grid = NetworkBuild.setup_grid(self.existing_grid)
            # merge the grid and nodes, returing the Net, Disjoint sets and RTree
            G, DS, R = NetworkBuild.grid_settlement_merge(grid, self._metric_graph())

        # Run the network optimization, filtering min_networks
        mst = nx.union_all(filter(lambda sub: len(sub.node) > 1,
                                  nx.connected_component_subgraphs(modBoruvka(G, DS, R))))
        nx.relabel_nodes(mst, {i: i+1 for i in mst.nodes()}, copy=False)

        # Add the existing grid to the dataset_store
        dataset_subnet = dataset_store.Subnet()
        for u, v in grid.edges():
            segment = dataset_store.Segment(u, v)
            segment.subnet_id = dataset_subnet.id
            segment.is_existing = True
            segment.weight = grid[u][v]['weight']
            self.dstore.session.add(segment)

        # Translate the NetworkX Graph to Dstore objects
        for subgraph in nx.connected_component_subgraphs(mst):
            # Initialize the subgraph in the datastore
            dataset_subnet = dataset_store.Subnet()
            self.dstore.session.add(dataset_subnet)
            self.dstore.session.commit()

            # Extend the dstore subnet with its segments
            for u, v, data in subgraph.edges(data=True):
                edge = u, v

                # If any fake nodes in the edge, add to the dstore
                for i, fake in enumerate([n for n in edge if mst.node[n]['mv'] == np.inf], 1):
                    dataset_node = self.dstore.addNode(subgraph.node[fake]['coords'], is_fake=True)
                    dataset_node.id = fake 
                    self.dstore.session.add(dataset_node)
                    self.dstore.session.commit()

                    # Edges should never be composed of two fake nodes
                    assert i <= 1

                # Add the edge to the subnet
                segment = dataset_store.Segment(*edge)
                segment.subnet_id = dataset_subnet.id
                segment.is_existing = False
                segment.weight = data['weight']
                self.dstore.session.add(segment)


        # Commit changes
        self.dstore.session.commit()
        self.mst = mst

    def _metric_graph(self):
        """Converts the dataset_store metrics to a nx graph"""

        data = [(i, {'mv': node.metric, 'coords': node.getCommonCoordinates()})
                 for i, node in enumerate(self.dstore.cycleNodes())]
        G = nx.Graph()
        G.add_nodes_from(data)
        return G

    def save_output(self):
        metric.saveMetricsConfigurationCSV(os.path.join(self.output, 'metrics-job-input'), self.metricConfiguration)
        metric.saveMetricsCSV(os.path.join(self.output, 'metrics-global'), self.metricModel, self.metricValueByOptionBySection)
        self.dstore.saveMetricsCSV(os.path.join(self.output, 'metrics-local'), self.metricModel, VS.HEADER_TYPE_ALIAS)
        self.dstore.saveSegmentsSHP(os.path.join(self.output, 'networks-proposed'), is_existing=False)


def dataset_store_to_nx_graph(dataset_store):
    
    all_nodes = list(dataset_store.cycleNodes()) + \
                list(dataset_store.cycleNodes(isFake=True))
    np_to_nx_id = {node.id: i for i, node in enumerate(all_nodes)} 

    data = [(i, {'mv': node.metric, 'coords': node.getCommonCoordinates(), 'np_id': node.id})
             for i, node in enumerate(all_nodes)]
    G = nx.Graph()
    G.add_nodes_from(data)

    seg_to_nx_ids = lambda seg:  (np_to_nx_id[seg.node1_id], np_to_nx_id[seg.node2_id])
    edges = [seg_to_nx_ids(s) for s in dataset_store.cycleSegments(is_existing=False)] 
    edge_weights = {seg_to_nx_ids(s): s.weight for s in dataset_store.cycleSegments(is_existing=False)} 
    edge_is_existing = {seg_to_nx_ids(s): s.is_existing for s in dataset_store.cycleSegments(is_existing=False)} 
    edge_subnet_id = {seg_to_nx_ids(s): s.subnet_id for s in dataset_store.cycleSegments(is_existing=False)} 
    G.add_edges_from(edges)
    nx.set_edge_attributes(G, 'weight', edge_weights)
    nx.set_edge_attributes(G, 'is_existing', edge_is_existing)
    nx.set_edge_attributes(G, 'subnet_id', edge_subnet_id)

    return G

if __name__ == '__main__':
    metrics = '/home/cjn/src/networkplanner/np/public/files/demographicsLL.csv'
    json_met = '/home/cjn/src/networkplanner/sample_metric_params.json'
    output_path = '/home/cjn/np_data/potou_500kwh_dmd_ModBo'
    grid = None # '/home/cjn/src/networkplanner/LeonaNetworksLL.shp' 

    netplan = NP_Boruvka_Interface(metrics, json_met, output_path, grid)
