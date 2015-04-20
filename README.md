# Networker

Networker contains tools for planning distribution networks from supply and demand data.  

It relies on numpy, networkx and also utilizes several spatial indexes.  

## Overview 

The main interfaces for generating networks are the following classes:
- NetworkerRunner:  Assumes a pre-computed budget for each node.  It follows 
networkx output conventions for spatial data. 
- NetworkPlannerRunner:  Uses existing networkplanner econometric models
to compute nodal budget (aka mvMax).  It's output follows networkplanner 
conventions.  

Included are wrapper scripts for leveraging these classes from the command line.  

## Examples

Once installed with source, navigate to the top level source directory and run:

```
> python scripts/run_networker.py test/networker_config_max100.json -o output
```

Once run, you should see an output directory with edges and nodes shapefiles 
(with supporting files).  This particular example had sufficient nodal "budget"
for the resulting network to be a minimum spanning tree.  

From within python, you can run this via the NetworkerRunner class and work
directly with the resulting GeoGraph (subtype of networkx Graph):

```
from networker import networker_runner
import json

# load demand nodes and generate geograph
cfg = json.load(open("test/networker_config_max100.json"))
nr = networker_runner.NetworkerRunner(cfg)
geograph = nr.build_network()

# plot it
import matplotlib.pyplot as plt
from networker import utils

utils.draw_geograph(geograph)
plt.show()
```

![geograph image](http://i.imgur.com/r7ei1VR.png)

If you want to leverage the econometric models from networkplanner to compute 
the budget (aka mvMax) values for nodes, use the NetworkPlannerRunner class.  
This will also generate outputs that are consistent with networkplanner.  
A script to run this is included.  Here's a sample run:

```
> python scripts/run_networkplanner.py test/networkplanner_config_leona_net.json -o output
```

Once complete, the output directory should contain a dataset.db and networks-proposed 
shapefiles (along with addtional networkplanner outputs).  

## Configuration

Both NetworkerRunner and NetworkPlannerRunner take json configuration files
which define their inputs and how they should run.  Annotated schemas for these
configuration files are:

[NetworkRunner Config](https://github.com/SEL-Columbia/networker/blob/master/networker/networker_config_schema.json)

[NetworkPlannerRunner Config](https://github.com/SEL-Columbia/networker/blob/master/networker/networkplanner_config_schema.json)

## Installation

The simplest method of installation is via [conda](http://www.continuum.io/blog/conda).  

Once conda is installed ([guide here](http://docs.continuum.io/anaconda/install.html)), 
setup a new python 2.7 environment via:

```
conda create -n networker python=2.7
```

Once created, `source activate networker` to activate the environment and 
install subsequent libraries in that environment.  

To use the NetworkPlannerRunner, you will need the networkplanner library 
installed in this environment:

```
# the geojson library will be installed from the ioos channel so add it
conda config --add channels 'ioos'
conda install -c sel networkplanner-metrics
```

Now install the networker library:

```
conda install -c sel networker
```

For running scripts and doing development, you will need the source locally so
just clone this repository.  

For development you will need to set your PYTHONPATH to the source location.  

## Testing

Once setup for development, run nosetests from the top level repo directory. 
