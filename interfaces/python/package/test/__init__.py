import unittest
from igraph.test import basic, structural, flow
from igraph import *

def suite():
    return unittest.TestSuite( \
        [basic.suite(),
	 structural.suite(),
	 flow.suite()] \
    )
    
def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
