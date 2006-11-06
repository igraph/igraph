import unittest
from igraph import *

class BasicTests(unittest.TestCase):
    def testGraphCreation(self):
	self.failUnless(isinstance(Graph(), Graph))
	
def suite():
    basic_suite = unittest.makeSuite(BasicTests)
    return unittest.TestSuite((basic_suite))

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

