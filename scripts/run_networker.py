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
        help="base directory which all paths are relative to")
parser.add_argument("--output_directory", "-o", \
        default=".", \
        help="directory where all output files will be written")

args = parser.parse_args()

cfg = json.load(open(args.config_filename))
# switch working dir AFTER loading config
os.chdir(args.working_directory)
nwk = networker_runner.NetworkRunner(cfg, args.output_directory)

try:
    nwk.validate()
except Exception as e:
    sys.exit("validation failed: {}".format(str(e)))

try:
    nwk.run()
except Exception as e:
    sys.exit("run failed: {}".format(str(e)))
