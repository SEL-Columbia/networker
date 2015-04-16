import argparse
import sys
import json
from np.lib import dataset_store
from networker import utils
from networker import networkplanner_runner


def get_geograph(filename):
    ds = dataset_store.load(filename)
    g = networkplanner_runner.dataset_store_to_geograph(ds)
    return g


parser = argparse.ArgumentParser(description=
        "Output networkplanner dataset store as networkx based js graph")
parser.add_argument("dataset_file", help="dataset file representing geograph")

args = parser.parse_args()

g = get_geograph(args.dataset_file)

js = utils.geograph_to_json(g)

json.dump(js, sys.stdout)
