#!/usr/bin/env python

import argparse
import json
import sys
import os
import logging
import networkx as nx
import networker
import networker.io as nio
from networker.classes.geograph import GeoGraph

# setup log
logger = logging.getLogger('networker')

logger.info("networker %s (Python %s)" % (
                networker.__version__,
                '.'.join(map(str, sys.version_info[:3]))))

parser = argparse.ArgumentParser(description="Project nodes to nearest point on network")
parser.add_argument("node_filename",
                    help="csv of nodes to be joined to network")
parser.add_argument("network_filename",
                    help="network nodes will be joined to")
parser.add_argument("--x_column", "-x", \
        default="x", \
        help="column name for x value in node csv")
parser.add_argument("--y_column", "-y", \
        default="y", \
        help="column name for x value in node csv")
parser.add_argument("--rtree", "-r", \
        dest="rtree", action="store_true",
        help="use rtree for segment lookup")
parser.add_argument("--output_directory", "-o", \
        default=".", \
        help="directory where all output files will be written")
parser.set_defaults(rtree=True)

args = parser.parse_args()

nodes = nio.load_nodes(args.node_filename,
                                    args.x_column,
                                    args.y_column)

net = nio.load_shp(args.network_filename, simplify=False)

# relabel nodes/coords so that merge can work
prefix = "net-"
nx.relabel_nodes(net,
    {n: prefix + str(n) for n in net.nodes()}, copy=False)
net.coords = {prefix + str(n): c for n, c in
    net.coords.items()}

def project_helper(use_rtree):
    if use_rtree:
        rtree = net.get_rtree_index()
        return net.project_onto(nodes, rtree_index=rtree)
    else:
        return net.project_onto(nodes)

# project_onto returns geograph with projected nodes PLUS
# the network edges they were projected onto
projected = project_helper(args.rtree)

# get only the projected edges
# construct geograph of only projected edges
projected_edges = GeoGraph(srs=projected.srs)
for node in nodes:
    edge = (node, projected.neighbors(node)[0])
    projected_edges.add_edge(*edge)
    projected_edges.coords[edge[0]] = projected.coords[edge[0]]
    projected_edges.coords[edge[1]] = projected.coords[edge[1]]

nio.write_shp(projected_edges, args.output_directory)

