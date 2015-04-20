NetworkBuild
============

An alternative algorithm for building an optimized infastructure network using nearest foreign neighbor queries.

Modified Borůvka's
--

This algorithm reduces the space complexity currently required in Network Planner by only considering edges between 
nearest foreign neighbors. This is made possible by a recursive space partitioning data structure offering cheap 
lookups in sublinear time. The core MST algorithm is largely influenced by Boruvka Meets Nearest Neighbors [1].


Spatial Index
--
Nodal coordinates are stored in a KDTree for efficient k-nearest neighbor lookups. Due to the periodicity in geodetic 
coordinate systems, the points are first projected into 3D cartesian space to preserve adjacency.

(example)[https://gfycat.com/RingedSarcasticFireant#]

[1] Mariano Tepper, Marta Mejail, Pablo Muse, Andres Almansa. Boruvka Meets Nearest Neighbors. 2011. <hal-00583120>

Motivation
==
Background
--
NetworkPlanner builds an optimized network that can be used in the initial design of extending electricity access in developing
nations. It accomplishes this by first computing a series of metrics, and then using those metrics to perform a contrained 
minumum spanning tree. These contraints include the two-way MV test, and the line (cycle) intersection test. Due to these
contraints however, the result is more likely a spanning forest. In practice this means the software uses a traditional minimum
spanning tree algorithm, namely [Kruskal's](http://en.wikipedia.org/wiki/Kruskal's_algorithm). 

Kruskal
--
```
Given G(V, E)
E = {u, v, dist(u, v) | u, v ϵ V}
Et = {}
for u, v in sorted(E, λ.u.v.d -> d):
    if subgraph[u] != subgraph[d]:
       subgraph[u] ⋃ subgraph[d]
       Et += (u, v)
MST(G) -> G(V, Et)
```
Kruskal's algorithm is extremely elegant and is easily extensible, such that the constraints NetworkPlanner requires can be
implemented with ease. You simply add the constraints beneath the requirement that the proposed edge spans disjoint components.
```
Given G(V, E)
E = {u, v, dist(u, v) | u, v ϵ V}
Et = {}
for u, v in sorted(E, λ.u.v.d -> d):
    if subgraph[u] != subgraph[d]:
      if successful two-way MV test:
         if u, v doesn't introduce cycles by intersection:
            subgraph[u] ⋃ subgraph[d]
            Et += (u, v)
NetOpt(G) -> G(V, Et)
```

This algorithm introduces a significant problem however when you want to optimize networks that include a large number of
nodes. The current NetworkPlanner will struggle on ```N ≅ 5000```. This is because the algorithm exhaustively checks every 
possible edge combination as evident in this line ```E = {u, v, dist(u, v) | u, v ϵ V}```. CS literature calls this ``O(n^2)`` 
space complexity, as the number of edges stored in memory scales with the square of the number of data points. Therefore with a 5000 node set, we are required to store 25,000,000 edges. Each of these edges has a corresponding weight that is used to 
prioritize their addition to the resulting MST. Each edge is likely a 64bit float, requiring ```(n^2*64) / 8  * 10^-9``` GB 
memory, or more concretely ~1.6GB for ```N=5000```. While many consumer PC's have this amount of memory, its apparent this is 
nearing the limits of whats possible.

Furthermore, this is just the baseline Kruskal algorithm not to mention the additional datastructures the software requires.
To further emphasize this limitation, consider the goal of optimizing 43,000 nodes, this would require ~120GB of memory!

Partitioning vs Optimization
---
The initial response to this problem is to simply chunk the problem into smaller sets and run the algorithm as is. In fact this
is what the team was manually doing for some time. This raises the problem that those we train to use NP are unable to 
replicate results that we present to them. On first instinct, we thought why not automate this partitioning? We talked through 
a fairly elegant solution to this, but we realized we would lose an aspect of the current algorithm. The Network propagation
effect which allows the spanning tree to make more expensive connections from a componenet as that component expands. By 
partitioning, we limit this effect. The optimal solution would be to reduce the ```O(n^2)``` space complexity by a factor of n 
such that we can run huge datasets. To do so we need to think about what distances we really need to compute. For instance 
```(n^2 / 2) - n``` edges can immediately be thrown out, as they are duplicates or zero ie edge ```(1, 2) == (2, 1)``` and edge
```(u, u) == (v, v) == 0```. Really, since we are interested in the minimal spanning tree we are most interested in the 
shortest edges or the edge between a node and its nearest neighbor.

(modified) Boruvka
--
[Boruvka](http://en.wikipedia.org/wiki/Borůvka's_algorithm) is the original MST algorithm, and with minor modifications fits 
our needs exactly. Boruvka works, by treating
everything as a disjoint component, and connecting it to its nearest disjoint neighbor. We can efficiently find these neighbors
using a KDTree with a modified traversal rule that allows queries restricted to subset of the data.
```
G(V, E)
DS = UnionFind()
Et, E' = [], PriorityQueue()
for v ∈ V:
  DS[v]
  dm, vm = FNN(v, V)
  push DS.queue[v] (v, vm) dm
while |Et| ≠ |V| - 1:
  for C ∈ components of DS:
    v, vm, dm = top DS.queue[C]
    while vm ∈ C:
      pop DS.queue[C]
      dm, um = FNN(v, {V - C})
      push DS.queue[v] (v, um) dm
      v, vm, dm = pop DS.queue[C]
    push E' (v, vm) dm
  while E':
    u, v, dm = pop E'
    if DS[u] ≠ DS[v]:
       DS[u] ⋃ DS[v]
       Et += (u, v)
MST(G) = G(V, Et)
```
This algorithm reduces the space complexity to ~O(n) as it only stores a reference from each node to its nearest foreign 
neighbor. Furthermore the runtime is reduced, as these priority queues are merged on each iteration of the outer loop and only 
the top of these queues are considered for addition to the outputted edge list. Distance computations are limited and only 
computed if the top of a queue is no longer a valid possible edge (edge connects to its own component, ie was optimal before 
merging components in the previous iteration). Updating the queues is managed in this inner loop:
```
while vm ∈ C:
      pop DS.queue[C]
      dm, um = FNN(v, {V - C})
      push DS.queue[v] (v, um) dm
      v, vm, dm = pop DS.queue[C]
    push E' (v, vm) dm
```
The factor that makes this algorithm possible and extensible to the NetworkPlanner constraints, is that like Kruskal its based 
on the premise of merging disjoint components and the logic can be implemented in much the same way. 
```
G(V, E)
DS = UnionFind()
Et, E' = [], PriorityQueue()
for v ∈ V:
  DS[v]
  dm, vm = FNN(v, V)
  push DS.queue[v] (v, vm) dm
while Et != state
  for C ∈ components of DS:
    v, vm, dm = top DS.queue[C]
    while vm ∈ C:
      pop DS.queue[C]
      dm, um = FNN(v, {V - C})
      push DS.queue[v] (v, um) dm
      v, vm, dm = pop DS.queue[C]
    push E' (v, vm) dm
  state = Et
  while E':
    u, v, dm = pop E'
    if DS[u] ≠ DS[v]:
      if successful two-way MV test:
         if u, v doesn't introduce cycles by intersection:
           DS[u] ⋃ DS[v]
           Et += (u, v)
NetOpt(G) = G(V, Et)
```
Other Optimizations and Changes
---
- Line Intersection Test
-- Shapely Intersection Test -> Efficient Linear Algebra Test
- Fake Node Projection
-- Shapely Line Interpolation -> Efficient Linear Algebra Projection
- MV updates
--Disjoint Grids cause global update in NP, these are considered seperate entities in NetworkBuild

Interfaces
==
NetworkBuild
--
The lightest interface to run this Network model is to take the metric output and existing network shapefile and use the 
NetworkBuild class.

```python
from networkbuild import NetworkBuild, modBoruvka

metric = csv
grid = shp
#reads in the data projects fake nodes and merges the points
interface = NetworkBuild(metric, Grid) 
#runs the network build algorithm
optnet = interface.build(modBoruvka, min_node_component=2) 
```
optnet is now a NetworkX object containing the optimized network

NetworkPlanner
---
If you have networkplanner on your machine you can run it with the new algorithm by using the included 
NetworkPlannerInterface.py, you can import this class and use it to run the metric model and then run the network build. This 
will result in an output matching the current NetworkPlanner behavior.
```python
from NetworkBuild import NP_Boruvka_Interface

input_nodes = '/Users/blogle/Desktop/demographicsXY/nodes.csv'
json_met = '/Users/blogle/Desktop/demographicsXY/metric_params.json'
grid = '/Users/blogle/Desktop/1979/networks-existing.shp'
output_path = '/Users/blogle/Desktop/test'

NP_Boruvka_Interface(input_nodes, json_met, grid, output_path)
```
