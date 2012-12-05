#!/usr/bin/env python
"""
Runs all the doctests in the igraph module
"""

import doctest
import igraph

if __name__ == "__main__":
    doctest.testmod(igraph)
    doctest.testmod(igraph.clustering)
    doctest.testmod(igraph.cut)
    doctest.testmod(igraph.datatypes)
    doctest.testmod(igraph.drawing.utils)
    doctest.testmod(igraph.formula)
    doctest.testmod(igraph.remote.nexus)
    doctest.testmod(igraph.statistics)
    doctest.testmod(igraph.utils)

