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
        self.failUnless(self.gempty.diameter(unconn=False) == 10)
        self.failUnless(self.g.diameter() == 2)
        self.failUnless(self.gdir.diameter(False) == 2)
        self.failUnless(self.gdir.diameter() == 3)
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
        self.failUnless(self.g.degree(loops=False) == [3, 3, 2, 2])
        self.failUnless(self.g.degree() == [5, 3, 2, 2])
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


class BiconnectedComponentTests(unittest.TestCase):
    g1 = Graph.Full(10)
    g2 = Graph(5, [(0,1),(1,2),(2,3),(3,4)])
    g3 = Graph(6, [(0,1),(1,2),(2,3),(3,0),(2,4),(2,5),(4,5)])

    def testBiconnectedComponents(self):
        s = self.g1.biconnected_components()
        self.failUnless(len(s) == 1 and s[0]==range(10))
        s, ap = self.g1.biconnected_components(True)
        self.failUnless(len(s) == 1 and s[0]==range(10))

        s = self.g3.biconnected_components()
        self.failUnless(len(s) == 2 and s[0]==[2,4,5] and s[1]==[0,1,2,3])
        s, ap = self.g3.biconnected_components(True)
        self.failUnless(len(s) == 2 and s[0]==[2,4,5] and \
          s[1]==[0,1,2,3] and ap == [2])

    def testArticulationPoints(self):
        self.failUnless(self.g1.articulation_points() == [])
        self.failUnless(self.g2.cut_vertices() == [1,2,3])
        self.failUnless(self.g3.articulation_points() == [2])


class CentralityTests(unittest.TestCase):
    def testEigenvectorCentrality(self):
        g = Graph.Star(11)
        cent = g.evcent()
        self.failUnless(cent.index(max(cent)) == 0)
        self.assertAlmostEquals(max(cent), 1.0, places=3)
        cent, ev = g.evcent(scale=False, return_eigenvalue=True)
        if cent[0]<0: cent = [-x for x in cent]
        self.failUnless(cent.index(max(cent)) == 0)
        self.assertAlmostEquals(cent[1]/cent[0], 0.3162, places=3)
        self.assertAlmostEquals(ev, 3.162, places=3)

class MiscTests(unittest.TestCase):
    def testConstraint(self):
        g = Graph(4, [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3)])
        self.failUnless(isinstance(g.constraint(), list)) # TODO check more

    def testTopologicalSorting(self):
        g = Graph(5, [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)], directed=True)
        self.failUnless(g.topological_sorting() == [0, 4, 1, 2, 3])
        self.failUnless(g.topological_sorting(IN) == [3, 4, 2, 1, 0])
        g.to_undirected()
        self.assertRaises(InternalError, g.topological_sorting)

    def testLineGraph(self):
        g = Graph(4, [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3)])
        el = g.linegraph().get_edgelist()
        el.sort()
        self.failUnless(el == [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (2, 4), (3, 4)])

        g = Graph(4, [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3)], directed=True)
        el = g.linegraph().get_edgelist()
        el.sort()
        self.failUnless(el == [(0, 2), (0, 4)])


def suite():
    simple_suite = unittest.makeSuite(SimplePropertiesTests)
    degree_suite = unittest.makeSuite(DegreeTests)
    local_transitivity_suite = unittest.makeSuite(LocalTransitivityTests)
    biconnected_suite = unittest.makeSuite(BiconnectedComponentTests)
    centrality_suite = unittest.makeSuite(CentralityTests)
    misc_suite = unittest.makeSuite(MiscTests)
    return unittest.TestSuite([simple_suite,
                               degree_suite,
                               local_transitivity_suite,
                               biconnected_suite,
                               centrality_suite,
                               misc_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

