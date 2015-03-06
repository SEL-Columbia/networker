# Comparitive Analysis

These are the results of an exhaustive comparitive analysis of the modKruskal 
based NetworkPlanner algorithm vs the modBoruvka based NetworkBuild algorithm.

The latter is being developed in effort to address performance issues and modularize the former.  

## Methodology

The results here are based on running NetworkPlanner and NetworkBuild on 
identical inputs and capturing performance statistics (run time, 
maximum memory, etc).  The outputs were then compared to determine 
correctness (assuming NetworkPlanner results are "correct").  

Scenarios were run with the following number of nodes:
100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200

These node sets were generated randomly following a uniform distribution of
latitude, longitude pairs between 0.0 and 1.0 degrees.  

Each scenario was then run with 2 demands:

- Medium (mvMax = median distance of all nearest neighbors)
- High (mvMax = max distance of all nearest neighbors)

These resulted in sparse (Medium) and more connected (High) networks.

Tests were performed on a rackspace hosted server with 15G of ram, 50G SSD 
drive and 2 vCPU's (the "15 GB Memory v1" flavor).  

## Results

The results are broken down into correctness and performance metrics.  


### Correctness

The following tables compare the results of running identical inputs
through NetworkPlanner (np) and NetworkBuild (nb).  

The results are *almost* identical.  Note that rows are missing where 
NetworkPlanner scenarios failed due to out of memory conditions.  

#### Medium Demand

| num_nodes | np_num_edges | nb_num_edges | np_sum_edge_lens   | nb_sum_edge_lens   |
|-----------|--------------|--------------|--------------------|--------------------|
| 100       | 32           | 32           | 119816.9153137207  | 119816.9153276966  |
| 200       | 64           | 64           | 175145.5360107422  | 175145.5354555251  |
| 400       | 123          | 123          | 225955.20516967773 | 225955.20591952585 |
| 800       | 258          | 258          | 337265.002243042   | 337265.00145902054 |
| 1600      | 509          | 509          | 449296.6186981201  | 449296.6191125113  |
| 3200      | 1050         | 1050         | 665813.1923933029  | 665813.1931482236  |
| 6400      | 2066         | 2066         | 945427.9955635071  | 945427.9959691941  | 


#### High Demand

| num_nodes | np_num_edges | nb_num_edges | np_sum_edge_lens   | nb_sum_edge_lens   |
|-----------|--------------|--------------|--------------------|--------------------|
| 100       | 98           | 98           | 722015.7960510254  | 722015.7931261827  |
| 200       | 198          | 198          | 1030542.8436889648 | 1030542.839258047  |
| 400       | 399          | 398          | 1423662.604522705  | 1421991.7194512442 |
| 800       | 799          | 797          | 2097236.5056762695 | 2088079.4905254904 |
| 1600      | 1599         | 1599         | 2914150.146613121  | 2914150.148561445  |
| 3200      | 3199         | 3190         | 4063716.528371811  | 4041834.9088378055 |


### Performance

The following graphs compare the performance of NetworkPlanner (np) and 
NetworkBuild (nb) in terms of run time and memory usage.  

The results contrast and show: 

- The quadratic memory usage of NetworkPlanner vs the linear memory usage of 
NetworkBuild.  

- The cubic(?) run time of NetworkPlanner vs the O(NlogN) run time of 
NetworkBuild.

- The fact that runs beyond 6400 nodes were not possible via NetworkPlanner due
to the high memory usage.  

Run times are in seconds, memory is in GB.  Blue represents the existing NetworkPlanner modKruskal based algorith.  Green represents the new NetworkBuild modBoruvka based algorithm.  

#### Medium Demand Run Time

![med_dmd_runtime](http://i.imgur.com/16Ausud.png)

#### High Demand Run Time

![high_dmd_runtime](http://i.imgur.com/SgZdbdf.png)

#### Medium Demand Memory 

![med_dmd_mem](http://i.imgur.com/TCeGe5H.png)

#### High Demand Memory 

![high_dmd_mem](http://i.imgur.com/cPg30tV.png)


