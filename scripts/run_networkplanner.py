#!/usr/bin/env python

import argparse
import json
import sys
import os
import networker
import logging
from networker import networkplanner_runner

# setup log
logger = logging.getLogger('networker')
logger.setLevel(logging.INFO)

ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

ch.setFormatter(formatter)

logging.getLogger('networker').propagate = False

logger.addHandler(ch)

logger.info("networker %s (Python %s)" % (
                networker.__version__,
                '.'.join(map(str, sys.version_info[:3]))))


parser = argparse.ArgumentParser(description="Run networkplanner")
parser.add_argument("config_filename",
                    help="json config file for running networkplanner")
parser.add_argument("--working_directory", "-w", \
        default=".", \
        help="base directory which all input paths are relative to")
parser.add_argument("--output_directory", "-o", \
        default=".", \
        help="directory where all output files will be written")

args = parser.parse_args()

# switch working dir BEFORE loading config
os.chdir(args.working_directory)
cfg = json.load(open(args.config_filename))
nwk = networkplanner_runner.NetworkPlannerRunner(cfg, args.output_directory)

try:
    nwk.validate()
except Exception as e:
    sys.exit("validation failed: {}".format(str(e)))

try:
    nwk.run()
except Exception as e:
    sys.exit("run failed: {}".format(str(e)))
