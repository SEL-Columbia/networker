#!/usr/bin/env python

import argparse
import sys
import os
import re
import logging
import functools
import networkx as nx
import networker
import networker.io as nio
from networker.exception import NetworkerException
from networker.classes.geograph import GeoGraph

# setup log
logger = logging.getLogger('networker')

logger.info("networker %s (Python %s)" % (
                networker.__version__,
                '.'.join(map(str, sys.version_info[:3]))))

parser = argparse.ArgumentParser(description="Compose GeoGraphs")
parser.add_argument("input_files", nargs="+",
                    help="files representing geographs to be composed/unioned (.shp, .csv, .json, .geojson)")
parser.add_argument("output", 
                    help="resulting compose destination (dir for .shp or .geojson)")
parser.add_argument("--x_column", "-x", \
                    default="X", \
                    help="column name for x value in node csv")
parser.add_argument("--y_column", "-y", \
                    default="Y", \
                    help="column name for x value in node csv")
parser.add_argument("--force_disjoint", "-f", action="store_true", default=False,
                    help="force nodes to be distinct across input files (by assigning new integer node ids) resulting in a disjoint union")
parser.add_argument("--match_radius", "-r", type=float, 
                    help="if specified, determines radius for merging nodes once compose has been performed (otherwise nodes are not merged based on proximity)")
args = parser.parse_args()

def union_reduce(left_geo, right_geo):
    return GeoGraph.compose(left_geo, right_geo, args.force_disjoint)

def read_wrap(filename):
    if re.search(r'\.csv$', filename):
        return nio.read_geograph(filename, args.x_column, args.y_column)
    else:
        return nio.read_geograph(filename)

geographs = []
for filename in args.input_files:
    geograph = read_wrap(filename)
    logger.info("filename %s, nodes: %s, edges: %s" % (filename, len(geograph.nodes()), len(geograph.edges())))
    geographs.append(geograph)

# geographs = (read_wrap(filename) for filename in args.input_files)

unioned_geograph = functools.reduce(union_reduce, geographs)
    
if args.match_radius is not None:
   unioned_geograph.merge_nearby_nodes(args.match_radius)

logger.info("output nodes: %s, edges: %s" % (len(unioned_geograph.nodes()), len(unioned_geograph.edges())))
nio.write_geograph(unioned_geograph, args.output)
