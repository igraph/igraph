import unittest
from igraph import *

class BipartiteTests(unittest.TestCase):
    def testCreateBipartite(self):
        g = Graph.Bipartite([0, 1]*5, [(0,1),(2,3),(4,5),(6,7),(8,9)])
        self.assertTrue(g.vcount() == 10 and g.ecount() == 5 and g.is_directed() == False)
        self.assertTrue(g.is_bipartite())
        self.assertTrue(g.vs["type"] == [False, True]*5)

    def testFullBipartite(self):
        g = Graph.Full_Bipartite(10, 5)
        self.assertTrue(g.vcount() == 15 and g.ecount() == 50 and g.is_directed() == False)
        expected = sorted([(i, j) for i in xrange(10) for j in xrange(10, 15)])
        self.assertTrue(sorted(g.get_edgelist()) == expected)
        self.assertTrue(g.vs["type"] == [False]*10 + [True]*5)

        g = Graph.Full_Bipartite(10, 5, directed=True, mode=OUT)
        self.assertTrue(g.vcount() == 15 and g.ecount() == 50 and g.is_directed() == True)
        self.assertTrue(sorted(g.get_edgelist()) == expected)
        self.assertTrue(g.vs["type"] == [False]*10 + [True]*5)

        g = Graph.Full_Bipartite(10, 5, directed=True, mode=IN)
        self.assertTrue(g.vcount() == 15 and g.ecount() == 50 and g.is_directed() == True)
        self.assertTrue(sorted(g.get_edgelist()) == sorted([(i,j) for j, i in expected]))
        self.assertTrue(g.vs["type"] == [False]*10 + [True]*5)

        g = Graph.Full_Bipartite(10, 5, directed=True)
        self.assertTrue(g.vcount() == 15 and g.ecount() == 100 and g.is_directed() == True)
        expected.extend([(j, i) for i, j in expected])
        expected.sort()
        self.assertTrue(sorted(g.get_edgelist()) == expected)
        self.assertTrue(g.vs["type"] == [False]*10 + [True]*5)

    def testIncidence(self):
        g = Graph.Incidence([[0, 1, 1], [1, 2, 0]])
        self.assertTrue(g.vcount() == 5 and g.ecount() == 4 and g.is_directed() == False)
        self.assertTrue(g.vs["type"] == [False]*2 + [True]*3)
        self.assertTrue(sorted(g.get_edgelist()) == [(0,3),(0,4),(1,2),(1,3)])

        g = Graph.Incidence([[0, 1, 1], [1, 2, 0]], multiple=True)
        self.assertTrue(g.vcount() == 5 and g.ecount() == 5 and g.is_directed() == False)
        self.assertTrue(g.vs["type"] == [False]*2 + [True]*3)
        self.assertTrue(sorted(g.get_edgelist()) == [(0,3),(0,4),(1,2),(1,3),(1,3)])

        g = Graph.Incidence([[0, 1, 1], [1, 2, 0]], directed=True)
        self.assertTrue(g.vcount() == 5 and g.ecount() == 4 and g.is_directed() == True)
        self.assertTrue(g.vs["type"] == [False]*2 + [True]*3)
        self.assertTrue(sorted(g.get_edgelist()) == [(0,3),(0,4),(1,2),(1,3)])

        g = Graph.Incidence([[0, 1, 1], [1, 2, 0]], directed=True, mode="in")
        self.assertTrue(g.vcount() == 5 and g.ecount() == 4 and g.is_directed() == True)
        self.assertTrue(g.vs["type"] == [False]*2 + [True]*3)
        self.assertTrue(sorted(g.get_edgelist()) == [(2,1),(3,0),(3,1),(4,0)])

    def testGetIncidence(self):
        mat = [[0, 1, 1], [1, 1, 0]]
        v1, v2 = [0, 1], [2, 3, 4]
        g = Graph.Incidence(mat)
        self.assertTrue(g.get_incidence() == (mat, v1, v2))
        g.vs["type2"] = g.vs["type"]
        self.assertTrue(g.get_incidence("type2") == (mat, v1, v2))
        self.assertTrue(g.get_incidence(g.vs["type2"]) == (mat, v1, v2))

    def testBipartiteProjection(self):
        g = Graph.Full_Bipartite(10, 5)

        g1, g2 = g.bipartite_projection()
        self.assertTrue(g1.isomorphic(Graph.Full(10)))
        self.assertTrue(g2.isomorphic(Graph.Full(5)))
        self.assertTrue(g.bipartite_projection(which=0).isomorphic(g1))
        self.assertTrue(g.bipartite_projection(which=1).isomorphic(g2))
        self.assertTrue(g.bipartite_projection(which=False).isomorphic(g1))
        self.assertTrue(g.bipartite_projection(which=True).isomorphic(g2))
        self.assertTrue(g1.es["weight"] == [5] * 45)
        self.assertTrue(g2.es["weight"] == [10] * 10)
        self.assertTrue(g.bipartite_projection_size() == (10, 45, 5, 10))

        g1, g2 = g.bipartite_projection(probe1=10)
        self.assertTrue(g1.isomorphic(Graph.Full(5)))
        self.assertTrue(g2.isomorphic(Graph.Full(10)))
        self.assertTrue(g.bipartite_projection(which=0).isomorphic(g2))
        self.assertTrue(g.bipartite_projection(which=1).isomorphic(g1))
        self.assertTrue(g.bipartite_projection(which=False).isomorphic(g2))
        self.assertTrue(g.bipartite_projection(which=True).isomorphic(g1))

        g1, g2 = g.bipartite_projection(multiplicity=False)
        self.assertTrue(g1.isomorphic(Graph.Full(10)))
        self.assertTrue(g2.isomorphic(Graph.Full(5)))
        self.assertTrue(g.bipartite_projection(which=0).isomorphic(g1))
        self.assertTrue(g.bipartite_projection(which=1).isomorphic(g2))
        self.assertTrue(g.bipartite_projection(which=False).isomorphic(g1))
        self.assertTrue(g.bipartite_projection(which=True).isomorphic(g2))
        self.assertTrue("weight" not in g1.edge_attributes())
        self.assertTrue("weight" not in g2.edge_attributes())

    def testIsBipartite(self):
        g = Graph.Star(10)
        self.assertTrue(g.is_bipartite() == True)
        self.assertTrue(g.is_bipartite(True) == (True, [False] + [True]*9))
        g = Graph.Tree(100, 3)
        self.assertTrue(g.is_bipartite() == True)
        g = Graph.Ring(9)
        self.assertTrue(g.is_bipartite() == False)
        self.assertTrue(g.is_bipartite(True) == (False, None))
        g = Graph.Ring(10)
        self.assertTrue(g.is_bipartite() == True)
        g += (2, 0)
        self.assertTrue(g.is_bipartite(True) == (False, None))
        
def suite():
    bipartite_suite = unittest.makeSuite(BipartiteTests)
    return unittest.TestSuite([bipartite_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

