import argparse
import numpy as np
import sys
from networkbuild.KDTree import KDTree
from networkbuild import utils

selectors = {
    "min": np.min, 
    "median": np.median, 
    "max": np.max 
}
                    
def generate_nodes(n, min_max, demand_selector="median"):
    assert n > 1, "Number of nodes must be > 1"

    coords = np.random.uniform(low=min_max[0], high=min_max[1], size=(n, 2))
    # proj coords into 3-space (center at center of earth)
    proj_coords = utils.cartesian_projection(coords)

    kdt = KDTree(proj_coords)

    # find all nearest neighbor distances (only use the index from the result tuple)
    index_set = set(range(len(coords)))
    nneigh_inds = [kdt.query_subset(proj_coords[i], \
                                    list(index_set - set([i])))[0] \
                   for i in range(len(coords))]
    nneigh_coords = coords[nneigh_inds]
    nn_dists = utils.get_hav_distance(coords[:, 0], coords[:, 1], 
                        nneigh_coords[:, 0], nneigh_coords[:, 1])

    # assign same demand value to all nodes based on nearest neighbor dists
    demand_vals = np.repeat(selectors[demand_selector](nn_dists), len(coords))
    coords_demand = np.hstack((coords, demand_vals[:, np.newaxis]))
    return coords_demand


parser = argparse.ArgumentParser(description=\
        "Generate a set of demand nodes with X, Y, Demand values.  \
         X, Y's are spherical lat/long, Demand is in meters")
parser.add_argument("num_nodes", type=int, help="number of nodes to generate")
parser.add_argument("--min_max", "-m", type=float, nargs=2, \
        help="range [MIN, MAX] of values of node x, y's", default=[0, 1])
parser.add_argument("--demand_selector", "-d", \
        choices=["min", "median", "max"], default="median", \
        help="method to select demand value from all nearest neighbor pairs")

args = parser.parse_args()

coords_demand = generate_nodes(args.num_nodes, args.min_max, args.demand_selector)
np.savetxt(sys.stdout, coords_demand, fmt='%.10f', delimiter=",")

