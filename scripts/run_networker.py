#!/usr/bin/env python

import argparse
import json
import sys
import os
from networker import networker_runner

parser = argparse.ArgumentParser(description="Run networker")
parser.add_argument("config_filename",
                    help="json config file for running networker")
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
nwk = networker_runner.NetworkerRunner(cfg, args.output_directory)

try:
    nwk.validate()
except Exception as e:
    sys.exit("validation failed: {}".format(str(e)))

try:
    nwk.run()
except Exception as e:
    sys.exit("run failed: {}".format(str(e)))
