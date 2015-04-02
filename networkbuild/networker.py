# -*- coding: utf-8 -*-

import networkx as nx
import pandas as pd

from networkbuild.utils import UnionFind
from networkbuild.classes import GeoGraph
from networkbuild.network_io import load_shp, write_shp

class Networker(object):

    """
    class for running minimum spanning forest algorithms on a set of
    spatially referenced nodes

    Attributes:
        config:  dict (potentially nested) of configuration params
            TODO:  detail the params here

    """

    def __init__(self, config):
        self.config = config

    def run(self):
        """
        run a minimum spanning forest algorithm on inputs and write output 
        based on configuration
        """



