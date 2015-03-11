import argparse
import sys
import json
from np.lib import dataset_store
from networkbuild import utils
from networkbuild import NetworkPlannerInterface as npi
                    
def get_nx_graph(filename):
    ds = dataset_store.load(filename)
    nxg = npi.dataset_store_to_nx_graph(ds)
    return nxg

parser = argparse.ArgumentParser(description=\
        "Output networkplanner dataset store as networkx based js graph")
parser.add_argument("dataset_file", help="dataset file representing nx graph")

args = parser.parse_args()

g = get_nx_graph(args.dataset_file)

js = utils.network_to_json(g)
 
json.dump(js, sys.stdout)

