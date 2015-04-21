# Modified Borůvka's

This algorithm reduces the space complexity of the legacy NetworkPlanner algorithm (Modified Kruskal's) by utilizing the concepts behind Boruvka's Minimum Spanning Tree algorithm.  This works was significantly influenced by Boruvka Meets Nearest Neighbors [1].

## Background
NetworkPlanner attempts to generate a least cost network via a constrained Minimum Spanning Tree algorithm. These contraints include the two-way MV test, and a line (cycle) intersection test. Due to these contraints, the result is more likely a spanning forest.  [Kruskal's Algorithm](http://en.wikipedia.org/wiki/Kruskal's_algorithm) was the original basis for this constrained Minimum Spanning algorithm.  

### Kruskal
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
Due to it's simplicity, Kruskal's algorithm makes it's extension to handle the required constraints relatively straigh-forward. 

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

This algorithm is O(n^2) memory complexity, making it impossible to run for large datasets (practically, those with > 10000 nodes).  


### (modified) Boruvka

[Boruvka](http://en.wikipedia.org/wiki/Borůvka's_algorithm) is the original MST algorithm, and with some modification (and some added complexity vs Kruskal's) fits 
our needs. Boruvka works by treating
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
### Other Optimizations and Changes

- Line Intersection Test
-- Shapely Intersection Test -> Efficient Linear Algebra Test
- Fake Node Projection
-- Shapely Line Interpolation -> Efficient Linear Algebra Projection
- MV updates
-- Disjoint Grids cause global update in NP, these are considered seperate entities in NetworkBuild

[1] Mariano Tepper, Marta Mejail, Pablo Muse, Andres Almansa. Boruvka Meets Nearest Neighbors. 2011. <hal-00583120>


