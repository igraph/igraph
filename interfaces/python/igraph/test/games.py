import unittest
from igraph import *

class GameTests(unittest.TestCase):
    def testGRG(self):
        g = Graph.GRG(50, 0.2)
        self.failUnless(isinstance(g, Graph))
        g = Graph.GRG(50, 0.2, True)
        self.failUnless(isinstance(g, Graph))
        self.failUnless("x" in g.vertex_attributes())
        self.failUnless("y" in g.vertex_attributes())
        self.failUnless(isinstance(Layout(zip(g.vs["x"], g.vs["y"])), Layout))

    def testForestFire(self):
        g=Graph.Forest_Fire(100, 0.1)
        self.failUnless(isinstance(g, Graph) and g.is_directed() == False)
        g=Graph.Forest_Fire(100, 0.1, directed=True)
        self.failUnless(isinstance(g, Graph) and g.is_directed() == True)

    def testRecentDegree(self):
        g=Graph.Recent_Degree(100, 5, 10)
        self.failUnless(isinstance(g, Graph))

    def testPreference(self):
        g=Graph.Preference(100, [1, 1], [[1, 0], [0, 1]])
        self.failUnless(isinstance(g, Graph) and len(g.clusters()) == 2)

        g=Graph.Preference(100, [1, 1], [[1, 0], [0, 1]], attribute="type")
        l=g.vs.get_attribute_values("type")
        self.failUnless(min(l) == 0 and max(l) == 1)

    def testAsymmetricPreference(self):
        g=Graph.Asymmetric_Preference(100, [[0, 1], [1, 0]], [[0, 1], [1, 0]])
        self.failUnless(isinstance(g, Graph) and len(g.clusters()) == 2)

        g=Graph.Asymmetric_Preference(100, [[0, 1], [1, 0]], [[1, 0], [0, 1]],\
                                      attribute="type")
        l=g.vs.get_attribute_values("type")
        l1=[i[0] for i in l]
        l2=[i[1] for i in l]
        self.failUnless(min(l1) == 0 and max(l1) == 1 and
                        min(l2) == 0 and max(l2) == 1)

        g=Graph.Asymmetric_Preference(100, [[0, 1], [1, 0]], [[1, 0], [0, 1]])
        self.failUnless(isinstance(g, Graph) and len(g.clusters()) == 1)

    def testWattsStrogatz(self):
        g=Graph.Watts_Strogatz(1, 20, 1, 0.2)
        self.failUnless(isinstance(g, Graph) and g.vcount()==20 and g.ecount()==20)

    def testRewire(self):
        # Undirected graph
        g=Graph.GRG(25, 0.4)
        degrees=g.degree()

        # Rewiring without loops
        g.rewire(10000)
        self.assertEquals(degrees, g.degree())
        self.assertTrue(g.is_simple())

        # Rewiring with loops (1)
        g.rewire(10000, mode="loops")
        self.assertEquals(degrees, g.degree())
        self.assertFalse(any(g.is_multiple()))

        # Rewiring with loops (2)
        g = Graph.Full(4)
        g[1,3] = 0
        degrees = g.degree()
        g.rewire(100, mode="loops")
        self.assertEquals(degrees, g.degree())
        self.assertFalse(any(g.is_multiple()))

        # Directed graph
        g=Graph.GRG(25, 0.4)
        g.to_directed("mutual")
        indeg, outdeg = g.indegree(), g.outdegree()
        g.rewire(10000)
        self.assertEquals(indeg, g.indegree())
        self.assertEquals(outdeg, g.outdegree())
        self.assertTrue(g.is_simple())

        # Directed graph with loops
        g.rewire(10000, mode="loops")
        self.assertEquals(indeg, g.indegree())
        self.assertEquals(outdeg, g.outdegree())
        self.assertFalse(any(g.is_multiple()))

def suite():
    game_suite = unittest.makeSuite(GameTests)
    return unittest.TestSuite([game_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

