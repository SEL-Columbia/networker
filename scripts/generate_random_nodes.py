import argparse
import numpy as np
import sys

def generate_nodes(n, x_range, y_range):
    assert n > 1, "Number of nodes must be > 1"

    x_coords = np.random.uniform(low=x_range[0], high=x_range[1], size=n)
    y_coords = np.random.uniform(low=y_range[0], high=y_range[1], size=n)
    return np.concatenate((x_coords[:, np.newaxis], y_coords[:, np.newaxis]),
                          axis=1)

parser = argparse.ArgumentParser(description=\
        "Generate a set of X, Y pairs uniformly distributed within "\
        "x_range and y_range")
parser.add_argument("num_nodes", type=int, help="number of nodes to generate")
parser.add_argument("--x_range", "-x", type=float, nargs=2, \
        help="range [MIN, MAX] of values of x coordinates", default=[0, 1])
parser.add_argument("--y_range", "-y", type=float, nargs=2, \
        help="range [MIN, MAX] of values of y coordinates", default=[0, 1])
args = parser.parse_args()

coords = generate_nodes(args.num_nodes, args.x_range, args.y_range)
np.savetxt(sys.stdout, coords, fmt='%.10f', delimiter=",")
