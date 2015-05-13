# Comparitive Analysis

Below is a comparitive analysis of the modKruskal algorithm vs the new 
modBoruvka algorithm for network planning.  

The latter has been developed to address performance issues and to better 
modularize the former.  

## Methodology

The results are based on running modKruskal and modBoruvka on 
identical inputs and capturing performance statistics (run time, 
maximum memory, etc).  The outputs were then compared to determine 
correctness (assuming the existing modKruskal results are "correct").  

Scenarios were run with the following numbers of nodes:
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

The following tables compare the metrics of running identical inputs
through modKruskal (mk) and modBoruvka (mb).  

The resulting networks are *almost* identical.  Note that modKruskal scenarios
failed with > 6400 Medium Demand nodes and > 3200 High Demand nodes due to the
large amount of memory it consumes.  

Edge lengths are in meters.

#### Medium Demand

| num_nodes | num_edges (mk) | num_edges (mb) | sum_edge_lens (mk) | sum_edge_lens (mb) |
|-----------|----------------|----------------|--------------------|--------------------|
| 100       | 32             | 32             | 119816             | 119816             |
| 200       | 64             | 64             | 175145             | 175145             |
| 400       | 123            | 123            | 225955             | 225955             |
| 800       | 258            | 258            | 337265             | 337265             |
| 1600      | 509            | 509            | 449296             | 449296             |
| 3200      | 1050           | 1050           | 665813             | 665813             |
| 6400      | 2066           | 2066           | 945427             | 945427             | 


#### High Demand

| num_nodes | num_edges (mk) | num_edges (mb) | sum_edge_lens (mk) | sum_edge_lens (mb) |
|-----------|----------------|----------------|--------------------|--------------------|
| 100       | 98             | 98             | 722015             | 722015             |
| 200       | 198            | 198            | 1030542            | 1030542            |
| 400       | 399 (1 unique) | 398            | 1423662            | 1421991            |
| 800       | 799 (2 unique) | 797            | 2097236            | 2088079            |
| 1600      | 1599           | 1599           | 2914150            | 2914150            |
| 3200      | 3199 (9 unique)| 3190           | 4063716            | 4041834            |


### Performance

The following graphs compare the performance of modKruskal (mk) and 
modBoruvka (mb) in terms of run time and memory usage.  

The results contrast and show: 

- The quadratic memory usage of modKruskal vs the linear memory usage of 
modBoruvka.  

- The cubic(?) run time of modKruskal vs the more linear (N*logN?) run time 
of modBoruvka.

- The fact that runs beyond 6400 nodes were not possible via modKruskal due
to the high memory usage.  

Run times are in seconds, memory is in GB.  Blue represents the existing 
modKruskal algorithm.  Green represents the new modBoruvka algorithm.  

#### Medium Demand Run Time

![med_dmd_runtime](http://i.imgur.com/h8DD01H.png)

#### High Demand Run Time

![high_dmd_runtime](http://i.imgur.com/3CHJN3V.png)

#### Medium Demand Memory 

![med_dmd_mem](http://i.imgur.com/32wrzBX.png)

#### High Demand Memory 

![high_dmd_mem](http://i.imgur.com/o2lnSnH.png)

The following are tabular performance comparisons of the Medium Demand and 
High demand scenarios.  

#### Medium Demand

| num_nodes | run_time (mk) | run_time (mb) | max_mem (mk) | mb_max_mem (mb) |
|-----------|---------------|---------------|--------------|-----------------|
| 100       | 3.08          | 4.53          | 70148        | 125908          |
| 200       | 3.44          | 7.67          | 76608        | 124504          |
| 400       | 7.12          | 14.22         | 92824        | 139892          |
| 800       | 17.26         | 29.96         | 120140       | 169392          |
| 1600      | 45.86         | 69.35         | 193300       | 231180          |
| 3200      | 136.58        | 155.48        | 572056       | 348300          |
| 6400      | 461.73        | 392.64        | 2059984      | 580980          |
| 12800     | NA            | 1193.24       | 7941136      | 1050076         |
| 25600     | NA            | 3447.14       | 15065732     | 2009748         |
| 51200     | NA            | 11126.55      | 15059432     | 3888820         |

#### High Demand

| num_nodes | run_time (mk) | run_time (mb) | max_mem (mk) | mb_max_mem (mb) |
|-----------|---------------|---------------|--------------|-----------------|
| 100       | 7.42          | 4.54          | 69284        | 125348          |
| 200       | 29.45         | 7.56          | 82728        | 124268          |
| 400       | 138.08        | 13.96         | 147536       | 140840          |
| 800       | 682.14        | 29.15         | 412824       | 172620          |
| 1600      | 3452.14       | 65.33         | 1467952      | 234876          |
| 3200      | 18428.88      | 150.69        | 5671536      | 344676          |
| 6400      | NA            | 363.22        | 15070520     | 588120          |
| 12800     | NA            | 953.82        | 15042516     | 1061492         |
| 25600     | NA            | 2695.52       | 15067320     | 2036224         |
| 51200     | NA            | 7569.01       | 15052820     | 3922480         |

