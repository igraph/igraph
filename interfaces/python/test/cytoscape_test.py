#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Simple script that tests CytoscapeGraphDrawer.

This script is kept separate from the unit tests as it is very
hard to test for the correctness of CytoscapeGraphDrawer without
a working instance of Cytoscape.

Prerequisites for running this test:

    1. Start Cytoscape
    2. Activate the Cytoscape RPC plugin, listening at port 9000

"""

from igraph import Graph
from igraph.drawing.graph import CytoscapeGraphDrawer

def test():
    g = Graph.GRG(100, 0.2)

    ### Adding network attributes
    g["name"] = "Network name"
    g["version"] = 5
    g["obsolete"] = False
    g["density"] = g.density()

    ### Adding vertex attributes
    # String attribute
    g.vs["name"] = ["Node %d" % (i+1) for i in xrange(g.vcount())]
    # Integer attribute
    g.vs["degree"] = g.degree()
    # Float attribute
    g.vs["pagerank"] = g.pagerank()
    # Boolean attribute
    g.vs["even"] = [i % 2 for i in xrange(g.vcount())]
    # Mixed attribute
    g.vs["mixed"] = ["abc", 123, None, 1.0] * ((g.vcount()+3) / 4)
    # Special attribute with Hungarian accents
    g.vs[0]["name"] = u"árvíztűrő tükörfúrógép ÁRVÍZTŰRŐ TÜKÖRFÚRÓGÉP"

    ### Adding edge attributes
    # String attribute
    g.es["name"] = ["Edge %d -- %d" % edge.tuple for edge in g.es]
    # Integer attribute
    g.es["multiplicity"] = g.count_multiple()
    # Float attribute
    g.es["betweenness"] = g.edge_betweenness()
    # Boolean attribute
    g.es["even"] = [i % 2 for i in xrange(g.ecount())]
    # Mixed attribute
    g.es["mixed"] = [u"yay", 123, None, 0.7] * ((g.ecount()+3) / 4)

    # Sending graph
    drawer = CytoscapeGraphDrawer()
    drawer.draw(g, layout="fr")

    # Fetching graph
    g2 = drawer.fetch()
    del g2.vs["hiddenLabel"]
    del g2.es["interaction"]

    # Check isomorphism
    result = g2.isomorphic(g)
    if not result:
        raise ValueError("g2 not isomorphic to g")

    # Check the graph attributes
    if set(g.attributes()) != set(g2.attributes()):
        raise ValueError("Graph attribute name set mismatch")
    for attr_name in g.attributes():
        if g[attr_name] != g2[attr_name]:
            raise ValueError("Graph attribute mismatch for %r" % attr_name)

    # Check the vertex attribute names
    if set(g.vertex_attributes()) != set(g2.vertex_attributes()):
        raise ValueError("Vertex attribute name set mismatch")

    # Check the edge attribute names
    if set(g.edge_attributes()) != set(g2.edge_attributes()):
        raise ValueError("Edge attribute name set mismatch")
        

if __name__ == "__main__":
    test()

