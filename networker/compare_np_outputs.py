import argparse
import networkx as nx
import sys
import pandas as pd
from np.lib import dataset_store
from networker import networkplanner


def get_nx_graph(filename):
    ds = dataset_store.load(filename)
    nxg = networkplanner.dataset_store_to_nx_graph(ds)
    return nxg


def num_edges(g):
    return len(g.edges())


def num_conn_comps(g):
    return len(list(nx.connected_components(g)))


def sum_edge_weights(g):
    weights = [e[2]['weight'] for e in g.edges(data=True)]
    return sum(weights)


def num_unconnected_nodes(g):
    return nx.degree_histogram(g)[0]


def edge_diff(g_left, g_right):
    get_edge_set = lambda g: set([frozenset(e) for e in g.edges()])
    edges_left = get_edge_set(g_left)
    edges_right = get_edge_set(g_right)
    only_in_left = [set(e) for e in edges_left if e not in edges_right]
    only_in_right = [set(e) for e in edges_right if e not in edges_left]
    return only_in_left, only_in_right


parser = argparse.ArgumentParser(description=
        "Compare 2 networkplanner outputs based on their dataset stores")
parser.add_argument("dataset_left", help="dataset in left most column")
parser.add_argument("dataset_right", help="dataset in right most column")

args = parser.parse_args()

g_left = get_nx_graph(args.dataset_left)
g_right = get_nx_graph(args.dataset_right)

edges_left, edges_right = edge_diff(g_left, g_right)

row_names = ['num_edges', 'num_unconnected_nodes', 'num_connected_components',
             'sum_edge_weights', 'num_uniq_edges', 'uniq_edges']

column_names = ['left', 'right']

df = pd.DataFrame(index=row_names, columns=column_names)

df.loc['num_edges']['left'] = num_edges(g_left)
df.loc['num_edges']['right'] = num_edges(g_right)
df.loc['num_unconnected_nodes']['left'] = num_unconnected_nodes(g_left)
df.loc['num_unconnected_nodes']['right'] = num_unconnected_nodes(g_right)
df.loc['num_connected_components']['left'] = num_conn_comps(g_left)
df.loc['num_connected_components']['right'] = num_conn_comps(g_right)
df.loc['sum_edge_weights']['left'] = sum_edge_weights(g_left)
df.loc['sum_edge_weights']['right'] = sum_edge_weights(g_right)
df.loc['num_uniq_edges']['left'] = len(edges_left)
df.loc['num_uniq_edges']['right'] = len(edges_right)
df.loc['uniq_edges']['left'] = edges_left
df.loc['uniq_edges']['right'] = edges_right

df.to_csv(sys.stdout)
