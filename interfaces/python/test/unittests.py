#!/usr/bin/env python
"""
Simple script that runs the unit tests of igraph.
"""

import sys
from igraph.test import run_tests

if __name__ == "__main__":
    if "-v" in sys.argv:
        verbosity = 2
    else:
        verbosity = 1
    run_tests(verbosity=verbosity)
