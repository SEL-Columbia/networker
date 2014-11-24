# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

#%pylab inline
import numpy as np
import networkx as nx
from networkbuild import NetworkBuild

from networkbuild import modBoruvka
from nose.tools import eq_

# <codecell>

def TestGrid():
    g_coords = np.array([[-7., 5],
                     [-5. , 0.],
                     [ 5. , 0.],
                     [7. , 5]])

    nodes = dict(enumerate(g_coords))
    grid = nx.Graph()
    grid.add_nodes_from(nodes)
    grid.add_edges_from(zip(grid.nodes()[:-1], grid.nodes()[1:]))
    
    nx.set_node_attributes(grid, 'coords', nodes)
    #prepend grid to labels
    grid = nx.relabel_nodes(grid, {n: 'grid-' + str(n) for n in grid.nodes()})
    # Set mv to 0
    nx.set_node_attributes(grid, 'mv', {n:0 for n in grid.nodes()})

    return grid.to_undirected()

def TestNet():
    coords = np.array([[0,  1], 
                       [0, 2.5],
                       [.5, 2.5],
                       [-.5, 2.5],
                       [0, -1]
                       ])
    
    mv = np.array([111196, # 0 grid connection
                   158000, 
                   60000,  # edges 1, 2, 3 create enough to connect to 0
                   60000,
                   300000  # 4 connects to grid providing subgraph[0] mv to connect 0 -> 1
                   ])
    
    net = nx.Graph()
    nodes = dict(enumerate(coords))
    net.add_nodes_from(nodes)
    nx.set_node_attributes(net, 'coords', nodes)
    nx.set_node_attributes(net, 'mv', dict(enumerate(mv)))
    
    return net.to_undirected()


def TestMSTBehavior():
    grid, net = TestGrid(), TestNet()

    #Note: Fake nodes integer label begins at the total number of nodes in net + 1
    #Hence why the fake node in the test is incremented by one on each iteration

    edges_at_iteration = [[(0, 1)], # 0 connects to fake_node
                          [(0, 2)], # 0, 1 cant connect 
                          [(0, 3), (1, 2)], # 1 connects to 2
                          [(4, 0), (1, 2), (1,3)], # 1 connects to 3
                          [(4, 5), (5, 0), (0, 1), 
                           (1, 2), (1, 3)]] #4 connects to grid providing mv for (1, 2)

    for n, _ in enumerate(net.node):
        sub = net.subgraph(range(n+1))
        G, DS, R = NetworkBuild.grid_settlement_merge(grid, sub)
        mst = modBoruvka(G, DS, R)
        #eq_(set(mst.edges()), set(edges_at_iteration[n]))

    """
    When further developing this test, it is helpful to graphically inspect the connections
    
    uncomment the below codeblock to plot the Network and Grid in ipython notebook with
    %pylab inline 
    magic enabled
    """
#     figsize(22, 6)
#     c = ['b' if x['mv'] == np.inf else 'r' for _, x in G.nodes(data=True)]
#     nx.draw_networkx(mst, nx.get_node_attributes(mst, 'coords'), node_color=c)
#     nx.draw_networkx(grid, nx.get_node_attributes(grid, 'coords'), node_color='m', edge_color='r')
