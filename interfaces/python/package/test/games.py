import unittest
from igraph import *

class GameTests(unittest.TestCase):
    def testGRG(self):
	g=Graph.GRG(50, 0.2)
	self.failUnless(isinstance(g, Graph))
	g=Graph.GRG(50, 0.2, True)
	self.failUnless(isinstance(g, Graph))

    def testRecentDegree(self):
	g=Graph.Recent_Degree(100, 5, 10)
	self.failUnless(isinstance(g, Graph))

def suite():
    game_suite = unittest.makeSuite(GameTests)
    return unittest.TestSuite((game_suite))

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

