import unittest
from igraph import *

class GeneratorTests(unittest.TestCase):
    def testStar(self):
        g=Graph.Star(5, "in")
        el=[(1,0),(2,0),(3,0),(4,0)]
        self.assertTrue(g.is_directed())
        self.assertTrue(g.get_edgelist() == el)
        g=Graph.Star(5, "out", center=2)
        el=[(2,0),(2,1),(2,3),(2,4)]
        self.assertTrue(g.is_directed())
        self.assertTrue(g.get_edgelist() == el)
        g=Graph.Star(5, "mutual", center=2)
        el=[(0,2),(1,2),(2,0),(2,1),(2,3),(2,4),(3,2),(4,2)]
        self.assertTrue(g.is_directed())
        self.assertTrue(sorted(g.get_edgelist()) == el)
        g=Graph.Star(5, center=3)
        el=[(0,3),(1,3),(2,3),(3,4)]
        self.assertTrue(not g.is_directed())
        self.assertTrue(sorted(g.get_edgelist()) == el)

    def testFamous(self):
        g=Graph.Famous("tutte")
        self.assertTrue(g.vcount() == 46 and g.ecount() == 69)
        self.assertRaises(InternalError, Graph.Famous, "unknown")

    def testFormula(self):
        tests = [
            (None, [], []),
            ("", [""], []),
            ("A", ["A"], []),
            ("A-B", ["A", "B"], [(0, 1)]),
            ("A --- B", ["A", "B"], [(0, 1)]),
            ("A--B, C--D, E--F, G--H, I, J, K",
                ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"],
                [(0,1), (2,3), (4,5), (6,7)]
            ),
            ("A:B:C:D -- A:B:C:D",
                ["A", "B", "C", "D"],
                [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)]
            ),
            ("A -> B -> C", ["A", "B", "C"], [(0,1), (1,2)]),
            ("A <- B -> C", ["A", "B", "C"], [(1,0), (1,2)]),
            ("A <- B -- C", ["A", "B", "C"], [(1,0)]),
            ("A <-> B <---> C <> D", ["A", "B", "C", "D"],
                [(0,1), (1,0), (1,2), (2,1), (2,3), (3,2)]),
            ("'this is' <- 'a silly' -> 'graph here'",
                ["this is", "a silly", "graph here"], [(1,0), (1,2)]),
            ("Alice-Bob-Cecil-Alice, Daniel-Cecil-Eugene, Cecil-Gordon",
                ["Alice", "Bob", "Cecil", "Daniel", "Eugene", "Gordon"],
                [(0,1),(1,2),(0,2),(2,3),(2,4),(2,5)]
            ),
            ("Alice-Bob:Cecil:Daniel, Cecil:Daniel-Eugene:Gordon",
                ["Alice", "Bob", "Cecil", "Daniel", "Eugene", "Gordon"],
                [(0,1),(0,2),(0,3),(2,4),(2,5),(3,4),(3,5)]
            ),
            ("Alice <-> Bob --> Cecil <-- Daniel, Eugene --> Gordon:Helen",
                ["Alice", "Bob", "Cecil", "Daniel", "Eugene", "Gordon", "Helen"],
                [(0,1),(1,0),(1,2),(3,2),(4,5),(4,6)]
            ),
            ("Alice -- Bob -- Daniel, Cecil:Gordon, Helen",
                ["Alice", "Bob", "Daniel", "Cecil", "Gordon", "Helen"],
                [(0,1),(1,2)]
            ),
            ('"+" -- "-", "*" -- "/", "%%" -- "%/%"',
                ["+", "-", "*", "/", "%%", "%/%"],
                [(0,1),(2,3),(4,5)]
            )
        ]
        for formula, names, edges in tests:
            g = Graph.Formula(formula)
            self.assertEqual(g.vs["name"], names)
            self.assertEqual(g.get_edgelist(), sorted(edges))

    def testFull(self):
        g=Graph.Full(20, directed=True)
        el=g.get_edgelist()
        el.sort()
        self.assertTrue(g.get_edgelist() == [(x, y) for x in range(20) for y in range(20) if x!=y])

    def testFullCitation(self):
        g=Graph.Full_Citation(20)
        el=g.get_edgelist()
        el.sort()
        self.assertTrue(not g.is_directed())
        self.assertTrue(el == [(x, y) for x in xrange(19) for y in xrange(x+1, 20)])

        g=Graph.Full_Citation(20, True)
        el=g.get_edgelist()
        el.sort()
        self.assertTrue(g.is_directed())
        self.assertTrue(el == [(x, y) for x in xrange(1, 20) for y in xrange(x)])

        self.assertRaises(InternalError, Graph.Full_Citation, -2)

    def testLCF(self):
        g1=Graph.LCF(12, (5, -5), 6)
        g2=Graph.Famous("Franklin")
        self.assertTrue(g1.isomorphic(g2))
        self.assertRaises(InternalError, Graph.LCF, 12, (5, -5), -3)

    def testKautz(self):
        g=Graph.Kautz(2, 2)
        deg_in=g.degree(mode=IN)
        deg_out=g.degree(mode=OUT)
        # This is not a proper test, but should spot most errors
        self.assertTrue(g.is_directed() and deg_in==[2]*12 and deg_out==[2]*12)

    def testDeBruijn(self):
        g=Graph.De_Bruijn(2, 3)
        deg_in=g.degree(mode=IN, loops=True)
        deg_out=g.degree(mode=OUT, loops=True)
        # This is not a proper test, but should spot most errors
        self.assertTrue(g.is_directed() and deg_in==[2]*8 and deg_out==[2]*8)

    def testSBM(self):
        pref_matrix = [[0.5, 0, 0], [0, 0, 0.5], [0, 0.5, 0]]
        n = 60
        types = [20, 20, 20]
        g = Graph.SBM(n, pref_matrix, types)

        # Simple smoke tests for the expected structure of the graph
        self.assertTrue(g.is_simple())
        self.assertFalse(g.is_directed())
        self.assertEqual([0]*20 + [1]*40, g.clusters().membership)
        g2 = g.subgraph(range(20, 60))
        self.assertTrue(not any(e.source // 20 == e.target // 20 for e in g2.es))

        # Check loops argument
        g = Graph.SBM(n, pref_matrix, types, loops=True)
        self.assertFalse(g.is_simple())
        self.assertTrue(sum(g.is_loop()) > 0)

        # Check directedness
        g = Graph.SBM(n, pref_matrix, types, directed=True)
        self.assertTrue(g.is_directed())
        self.assertTrue(sum(g.is_mutual()) < g.ecount())
        self.assertTrue(sum(g.is_loop()) == 0)

        # Check error conditions
        self.assertRaises(InternalError, Graph.SBM, -1, pref_matrix, types)
        self.assertRaises(InternalError, Graph.SBM, 61, pref_matrix, types)
        pref_matrix[0][1] = 0.7
        self.assertRaises(InternalError, Graph.SBM, 60, pref_matrix, types)

    def testWeightedAdjacency(self):
        mat = [[0, 1, 2, 0], [2, 0, 0, 0], [0, 0, 2.5, 0], [0, 1, 0, 0]]

        g = Graph.Weighted_Adjacency(mat, attr="w0")
        el = g.get_edgelist()
        self.assertTrue(el == [(0,1), (0,2), (1,0), (2,2), (3,1)])
        self.assertTrue(g.es["w0"] == [1, 2, 2, 2.5, 1])

        g = Graph.Weighted_Adjacency(mat, mode="plus")
        el = g.get_edgelist()
        self.assertTrue(el == [(0,1), (0,2), (1,3), (2,2)])
        self.assertTrue(g.es["weight"] == [3, 2, 1, 2.5])

        g = Graph.Weighted_Adjacency(mat, attr="w0", loops=False)
        el = g.get_edgelist()
        self.assertTrue(el == [(0,1), (0,2), (1,0), (3,1)])
        self.assertTrue(g.es["w0"] == [1, 2, 2, 1])

        
def suite():
    generator_suite = unittest.makeSuite(GeneratorTests)
    return unittest.TestSuite([generator_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

