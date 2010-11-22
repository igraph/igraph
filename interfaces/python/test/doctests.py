#!/usr/bin/env python
"""
Runs all the doctests in the igraph module
"""

import doctest
import igraph

if __name__ == "__main__":
    doctest.testmod(igraph)
    doctest.testmod(igraph.statistics)

