# vim:ts=4 sw=4 sts=4:
import unittest
from igraph import *

class GraphAdjacencyMatrixLikeIndexingTests(unittest.TestCase):
    def testSingleEdgeRetrieval(self):
        g = Graph.Famous("krackhardt_kite")
        for v1, v2 in g.get_edgelist():
            self.assertEqual(g[v1, v2], 1)
            self.assertEqual(g[v2, v1], 1)
        for v1 in xrange(g.vcount()):
            for v2 in set(range(g.vcount())) - set(g.neighbors(v1)):
                self.assertEqual(g[v1, v2], 0)
                self.assertEqual(g[v2, v1], 0)

        g.add_edge(1, 1)
        self.assertEqual(g[1, 1], 1)

    def testSingleEdgeRetrievalWeights(self):
        g = Graph.Famous("krackhardt_kite")
        g.es["weight"] = range(g.ecount())
        for idx, (v1, v2) in enumerate(g.get_edgelist()):
            self.assertEqual(g[v1, v2], idx)
            self.assertEqual(g[v2, v1], idx)
        for v1 in xrange(g.vcount()):
            for v2 in set(range(g.vcount())) - set(g.neighbors(v1)):
                self.assertEqual(g[v1, v2], 0)
                self.assertEqual(g[v2, v1], 0)

    def testSingleEdgeRetrievalAttrName(self):
        g = Graph.Famous("krackhardt_kite")
        g.es["value"] = range(20, g.ecount()+20)
        for idx, (v1, v2) in enumerate(g.get_edgelist()):
            self.assertEqual(g[v1, v2, "value"], idx+20)
            self.assertEqual(g[v2, v1, "value"], idx+20)
        for v1 in xrange(g.vcount()):
            for v2 in set(range(g.vcount())) - set(g.neighbors(v1)):
                self.assertEqual(g[v1, v2, "value"], 0)
                self.assertEqual(g[v2, v1, "value"], 0)


def suite():
    adjacency_suite = unittest.makeSuite(GraphAdjacencyMatrixLikeIndexingTests)
    return unittest.TestSuite([adjacency_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()
