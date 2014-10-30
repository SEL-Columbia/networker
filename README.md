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
Nodal coordinates are stored in a KDTree for efficient k-nearest neighbor lookups. Due to the periodicity in geodetic coordinate systems, the points are first projected into 3D cartesian space to preserve adjacency.

![https://gfycat.com/RingedSarcasticFireant#](http://giant.gfycat.com/RingedSarcasticFireant.gif)

[1] Mariano Tepper, Marta Mejail, Pablo Muse, Andres Almansa. Boruvka Meets Nearest Neighbors. 2011. <hal-00583120>
