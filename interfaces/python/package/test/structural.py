import unittest
from igraph import *


class SimplePropertiesTests(unittest.TestCase):
    gfull  = Graph.Full(10)
    gempty = Graph(10)
    g = Graph(4, [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3)])
    gdir = Graph(4, [(0, 1), (0, 2), (1, 2), (2, 1), (0, 3), (1, 3), (3, 0)], directed=True)
    tree = Graph.Tree(10, 3)

    def testDensity(self):
	self.failUnless(self.gfull.density() == 1.0)
	self.failUnless(self.gempty.density() == 0.0)
	self.failUnless(self.g.density() == 5.0/6.0)
	self.failUnless(self.g.density(True) == 5.0/8.0)
	self.failUnless(self.gdir.density() == 7.0/12.0)
	self.failUnless(self.gdir.density(True) == 7.0/16.0)
	self.failUnless(self.tree.density() == 0.2)
	
    def testDiameter(self):
	self.failUnless(self.gfull.diameter() == 1)
	self.failUnless(self.gempty.diameter() == 10)
	self.failUnless(self.g.diameter() == 2)
	self.failUnless(self.gdir.diameter() == 2)
	self.failUnless(self.gdir.diameter(True) == 3)
	self.failUnless(self.tree.diameter() == 4)
	
    def testTransitivity(self):
	self.failUnless(self.gfull.transitivity_undirected() == 1.0)
	self.failUnless(self.tree.transitivity_undirected() == 0.0)
	self.failUnless(self.g.transitivity_undirected() == 0.75)


class DegreeTests(unittest.TestCase):
    gfull  = Graph.Full(10)
    gempty = Graph(10)
    g = Graph(4, [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (0, 0)])
    gdir = Graph(4, [(0, 1), (0, 2), (1, 2), (2, 1), (0, 3), (1, 3), (3, 0)], directed=True)
    tree = Graph.Tree(10, 3)

    def testDegree(self):
	self.failUnless(self.gfull.degree() == [9] * 10)
	self.failUnless(self.gempty.degree() == [0] * 10)
	self.failUnless(self.g.degree() == [3, 3, 2, 2])
	self.failUnless(self.g.degree(loops=True) == [5, 3, 2, 2])
	self.failUnless(self.gdir.degree(type=IN) == [1, 2, 2, 2])
	self.failUnless(self.gdir.degree(type=OUT) == [3, 2, 1, 1])
	self.failUnless(self.gdir.degree(type=ALL) == [4, 4, 3, 3])

    def testMaxDegree(self):
	self.failUnless(self.gfull.maxdegree() == 9)
	self.failUnless(self.gempty.maxdegree() == 0)
	self.failUnless(self.g.maxdegree() == 3)
	self.failUnless(self.g.maxdegree(loops=True) == 5)
	self.failUnless(self.g.maxdegree([1, 2], loops=True) == 3)
	self.failUnless(self.gdir.maxdegree(type=IN) == 2)
	self.failUnless(self.gdir.maxdegree(type=OUT) == 3)
	self.failUnless(self.gdir.maxdegree(type=ALL) == 4)
	
	

class LocalTransitivityTests(unittest.TestCase):
    def testLocalTransitivityFull(self):
	trans = Graph.Full(10).transitivity_local_undirected()
	self.failUnless(trans == [1.0]*10)
	
    def testLocalTransitivityTree(self):
	trans = Graph.Tree(10, 3).transitivity_local_undirected()
	self.failUnless(trans[0:3] == [0.0, 0.0, 0.0])

    def testLocalTransitivityHalf(self):
	g = Graph(4, [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3)])
	trans = g.transitivity_local_undirected()
	trans = [round(x, 3) for x in trans]
	self.failUnless(trans == [0.667, 0.667, 1.0, 1.0])

    def testLocalTransitivityPartial(self):
	g = Graph(4, [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3)])
	trans = g.transitivity_local_undirected([1,2])
	trans = [round(x, 3) for x in trans]
	self.failUnless(trans == [0.667, 1.0])

def suite():
    simple_suite = unittest.makeSuite(SimplePropertiesTests)
    degree_suite = unittest.makeSuite(DegreeTests)
    local_transitivity_suite = unittest.makeSuite(LocalTransitivityTests)
    return unittest.TestSuite((simple_suite,
			       degree_suite,
	                       local_transitivity_suite))

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

