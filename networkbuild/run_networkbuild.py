import argparse
from networkbuild.NetworkPlannerInterface import NP_Boruvka_Interface

parser = argparse.ArgumentParser(description="Run NetworkBuilder")
parser.add_argument("demand_file_name",
                    help="demand input file with X, Y and 'metric' column")
parser.add_argument('metric_model_params',
                    help="model parameters json file")
parser.add_argument('output_path', help="directory where outputs will be placed")
parser.add_argument("--network_file_name", "-n",
                    help="existing network input file")
                   
args = parser.parse_args()

# Run NP Boruvka and output results in appropriate dir
NP_Boruvka_Interface(args.demand_file_name, args.metric_model_params, \
                     args.output_path, args.network_file_name)
