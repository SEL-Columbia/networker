import argparse
import sys
import numpy as np
import networkx as nx
import networker.io as nio
from networker.classes.geograph import GeoGraph
import networker.geomath as gm


def generate_mst(coords):
    """ 
    Generate a min spanning tree based on coordinate 
    distances
    """
    input_proj = gm.PROJ4_LATLONG
    if gm.is_in_lon_lat(coords):
        input_proj = gm.PROJ4_LATLONG
    else:
        input_proj = gm.PROJ4_FLAT_EARTH

    node_dict = dict(enumerate(coords))

    geo_nodes = GeoGraph(input_proj, node_dict)
    geo_full = geo_nodes.get_connected_weighted_graph()
    geo_mst = nx.minimum_spanning_tree(geo_full)
    geo_nodes.add_edges_from(geo_mst.edges(data=True))
    return geo_nodes


parser = argparse.ArgumentParser(description=\
        "Generate minimum spanning tree from input nodes\n"\
        "takes x,y's (no column header) from stdin\n"\
        "writes json graph to stdout")

args = parser.parse_args()

nodes = np.loadtxt(sys.stdin, delimiter=",", ndmin=2)
geograph = generate_mst(nodes)
nio.write_json(geograph, sys.stdout)
