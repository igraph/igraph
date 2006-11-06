import unittest
from igraph.test import basic
from igraph import *

def suite():
    return unittest.TestSuite( \
        basic.suite() \
    )
    
def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
