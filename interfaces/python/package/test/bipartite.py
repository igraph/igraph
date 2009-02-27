import unittest
from igraph import *

class BipartiteTests(unittest.TestCase):
    def testIsBipartite(self):
        g = Graph.Star(10)
        self.failUnless(g.is_bipartite() == True)
        self.failUnless(g.is_bipartite(True) == (True, [False] + [True]*9))
        g = Graph.Tree(100, 3)
        self.failUnless(g.is_bipartite() == True)
        g = Graph.Ring(9)
        self.failUnless(g.is_bipartite() == False)
        self.failUnless(g.is_bipartite(True) == (False, None))
        g = Graph.Ring(10)
        self.failUnless(g.is_bipartite() == True)
        g += (2, 0)
        self.failUnless(g.is_bipartite(True) == (False, None))
        
def suite():
    bipartite_suite = unittest.makeSuite(BipartiteTests)
    return unittest.TestSuite([bipartite_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

