import unittest
from igraph import *

class CliqueTests(unittest.TestCase):
    def setUp(self):
	self.g=Graph.Full(6)
	self.g.delete_edges([(0, 1), (0, 2), (3, 5)])

    def testCliques(self):
	tests = {(4, -1): [[1, 2, 3, 4], [1, 2, 4, 5]],
		 (2, 2): [[0, 3], [0, 4], [0, 5],
			  [1, 2], [1, 3], [1, 4], [1, 5],
			  [2, 3], [2, 4], [2, 5], [3, 4], [4, 5]],
		 (-1, -1): [[0], [1], [2], [3], [4], [5],
			    [0, 3], [0, 4], [0, 5],
			    [1, 2], [1, 3], [1, 4], [1, 5],
			    [2, 3], [2, 4], [2, 5], [3, 4], [4, 5],
			    [0, 3, 4], [0, 4, 5],
			    [1, 2, 3], [1, 2, 4], [1, 2, 5],
			    [1, 3, 4], [1, 4, 5], [2, 3, 4], [2, 4, 5],
			    [1, 2, 3, 4], [1, 2, 4, 5]]}
	for (lo, hi), exp in tests.iteritems():
	    self.assertEqual(map(tuple, exp), self.g.cliques(lo, hi))

    def testLargestCliques(self):
	self.assertEqual(self.g.largest_cliques(),
			 [(1, 2, 3, 4), (1, 2, 4, 5)])

    def testMaximalCliques(self):
	self.assertEqual(self.g.maximal_cliques(),
			 [(0, 3, 4), (0, 4, 5),
			  (1, 2, 3, 4), (1, 2, 4, 5)])

    def testCliqueNumber(self):
	self.assertEqual(self.g.clique_number(), 4)
	self.assertEqual(self.g.omega(), 4)

class IndependentVertexSetTests(unittest.TestCase):
    def setUp(self):
	self.g1=Graph.Tree(5, 2, TREE_OUT)
	self.g2=Graph.Tree(10, 2, TREE_OUT)

    def testIndependentVertexSets(self):
	tests = {(4, -1): [],
		 (2, 2): [(0, 3), (0, 4), (1, 2), (2, 3), (2, 4), (3, 4)],
		 (-1, -1): [(0,), (1,), (2,), (3,), (4,),
			    (0, 3), (0, 4), (1, 2), (2, 3), (2, 4),
			    (3, 4), (0, 3, 4), (2, 3, 4)]}
	for (lo, hi), exp in tests.iteritems():
	    self.assertEqual(exp, self.g1.independent_vertex_sets(lo, hi))

    def testLargestIndependentVertexSets(self):
	self.assertEqual(self.g1.largest_independent_vertex_sets(),
			 [(0, 3, 4), (2, 3, 4)])

    def testMaximalIndependentVertexSets(self):
	self.assertEqual(self.g2.maximal_independent_vertex_sets(),
			 [(0, 3, 4, 5, 6), (0, 3, 5, 6, 9),
			  (0, 4, 5, 6, 7, 8), (0, 5, 6, 7, 8, 9),
			  (1, 2, 7, 8, 9), (1, 5, 6, 7, 8, 9),
			  (2, 3, 4), (2, 3, 9), (2, 4, 7, 8)])

    def testIndependenceNumber(self):
	self.assertEqual(self.g2.independence_number(), 6)
	self.assertEqual(self.g2.alpha(), 6)

def suite():
    clique_suite = unittest.makeSuite(CliqueTests)
    indvset_suite = unittest.makeSuite(IndependentVertexSetTests)
    return unittest.TestSuite([clique_suite, indvset_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

