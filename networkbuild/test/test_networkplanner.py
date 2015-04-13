# -*- coding: utf-8 -*-

import os, itertools, json
import numpy as np
import networkx as nx
import networkbuild.classes as cls
from networkbuild import networkplanner
from networkbuild import network_io
from networkbuild import geo_math as gm 
from networkbuild.algorithms import mod_boruvka

from nose.tools import eq_


def networkplanner_run(config_file, known_results_file):
    # get config and run
    cfg_path = os.path.join(os.path.dirname(\
        os.path.abspath(__file__)), config_file)
    cfg = json.load(open(cfg_path))
    nwk_p = networkplanner.NetworkPlanner(cfg)
    nwk_p.run()

    # compare this run against existing results
    test_geo = network_io.load_shp(os.path.join(cfg['output_directory'], \
        "networks-proposed.shp"))
    known_geo = network_io.load_shp(known_results_file)
    # compare sets of edges

    test_edges = test_geo.get_coord_edge_set()
    known_edges = known_geo.get_coord_edge_set()

    assert test_edges == known_edges, \
        "edges in test do not match known results"

    # This is redundant (but leave it here as an example as it's more lenient
    # than test above)
    # assert nx.is_isomorphic(test_geo, known_geo), \
    #    "test result graph is not isomorphic to known result graph"

def test_networkplanner_run():
    """ test on randomly generated set of nodes (demand only) """

    run_config = "networkplanner_config_pop100.json"
    results_file = "data/pop_100/networks-proposed.shp"

    networkplanner_run(run_config, results_file)

def test_networkplanner_leona_run():
    """ test on randomly generated set of nodes (demand only) """

    run_config = "networkplanner_config_leona_net.json"
    results_file = "data/leona/expected/networks-proposed.shp"

    networkplanner_run(run_config, results_file)
 
