"""
Run through of simple networkplanner scenario
"""
import numpy as np
import pandas as pd
import networkx as nx
from networker.classes.geograph import GeoGraph
from networker import utils, networker_runner

"""
Given a set of \"Demand Nodes\" over a region

For simplicity, assume all have 0 supply
"""

def random_points(num_points, upto):
    x = np.random.randint(upto, size=num_points)
    y = np.random.randint(upto, size=num_points)
    return np.array([x, y]).T

NUM_NODES = 10
REGION_WIDTH = 1000

coords = random_points(NUM_NODES, REGION_WIDTH)
# TODO:  vary the distribution
households = np.random.randint(20, high=10000, size=10)
demand_nodes = pd.DataFrame({'households': households})

# TODO:  plot them


# Assign costs
"""
## Simple Generation and Distribution costs 

for $node_k$:

$$generation_k = households_k * 1500$$
$$distribution_k = households_k * 500$$
$$standalone_k = generation_k + distribution_k$$

$$total = \sum_{k=0}^n standalone_k$$
"""

GEN_COST_PER_HH = 1500
DISTR_COST_PER_HH = 500

demand_nodes.generation = demand_nodes.households * GEN_COST_PER_HH
demand_nodes.distribution = demand_nodes.households * DISTR_COST_PER_HH
total_generation = sum(demand_nodes.generation)
average_generation = total_generation / NUM_NODES
total_standalone = sum(demand_nodes.generation + demand_nodes.distribution)
average_standalone = total_standalone / NUM_NODES
print("total standalone: %s" % total_standalone)

# TODO:  plot standalone costs (map plus distribution)

"""
## What if we share via a network

Cheapest way to connect everyone is MST
"""

geo_nodes = GeoGraph(coords=coords)
geo_full = geo_nodes.get_connected_weighted_graph()
geo_mst = nx.minimum_spanning_tree(geo_full)
geo_nodes.add_edges_from(geo_mst.edges(data=True))
# utils.draw_geograph(geo_nodes)

# set the line_cost_per_km so that standalone is cheaper for some
total_line_length = sum(e[2]['weight'] for e in geo_mst.edges(data=True))
mean_inter_node_dist = total_line_length / NUM_NODES
line_cost_per_km = average_generation / mean_inter_node_dist
total_grid = total_line_length * line_cost_per_km
print("total grid: %s" % total_grid)

"""

Assume 
$$grid_k = linelength_k * unitlinecost$$

Now
$$total = \sum_{k=0}^n standalone_k[standalone_k < grid_k] + grid_k[grid_k < standalone_k]$$

We want to minimize this
"""

"""
## Our approach "modMST"
"""

demand_geograph = GeoGraph(coords=coords)
for node in demand_geograph.nodes_iter():
    demand_geograph.node[node] = {'budget': demand_nodes.generation[node] / line_cost_per_km}

msf = networker_runner.build_network(demand_geograph)
# utils.draw_geograph(msf)

total_line_length = sum(e[2]['weight'] for e in msf.edges(data=True))
total_mod_grid = total_line_length * line_cost_per_km

non_grid_nodes = [node for node in range(NUM_NODES) if node not in msf.edge]
total_mod_standalone = sum(demand_nodes.generation[node] + 
                           demand_nodes.distribution[node] 
                           for node in non_grid_nodes)

total_mod = total_mod_standalone + total_mod_grid
print("standalone: %s, grid %s, total: %s" % (total_mod_standalone, total_mod_grid, total_mod))
