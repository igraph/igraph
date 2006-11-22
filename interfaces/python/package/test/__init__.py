import unittest
from igraph.test import basic, structural
from igraph import *

def suite():
    return unittest.TestSuite( \
        [basic.suite(),
	structural.suite()] \
    )
    
def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
