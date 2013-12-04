import unittest
from igraph import *

class OperatorTests(unittest.TestCase):
    def testMultiplication(self):
        g = Graph.Full(3)*3
        self.assertTrue(g.vcount() == 9 and g.ecount() == 9
                        and g.clusters().membership == [0,0,0,1,1,1,2,2,2])

    def testIntersection(self):
        g = Graph.Tree(7, 2) & Graph.Lattice([7])
        self.assertTrue(g.get_edgelist() == [(0, 1)])

    def testUnion(self):
        g = Graph.Tree(7, 2) | Graph.Lattice([7])
        self.assertTrue(g.vcount() == 7 and g.ecount() == 12)

    def testInPlaceAddition(self):
        g = Graph.Full(3)
        orig = g

        # Adding vertices
        g += 2
        self.assertTrue(g.vcount() == 5 and g.ecount() == 3
                        and g.clusters().membership == [0,0,0,1,2])

        # Adding a vertex by name
        g += "spam"
        self.assertTrue(g.vcount() == 6 and g.ecount() == 3
                        and g.clusters().membership == [0,0,0,1,2,3])

        # Adding a single edge
        g += (2, 3)
        self.assertTrue(g.vcount() == 6 and g.ecount() == 4
                        and g.clusters().membership == [0,0,0,0,1,2])

        # Adding two edges
        g += [(3, 4), (2, 4), (4, 5)]
        self.assertTrue(g.vcount() == 6 and g.ecount() == 7
                        and g.clusters().membership == [0]*6)

        # Adding two more vertices
        g += ["eggs", "bacon"]
        self.assertEqual(g.vs["name"], [None, None, None, None, None,
            "spam", "eggs", "bacon"])

        # Did we really use the original graph so far?
        # TODO: disjoint union should be modified so that this assertion
        # could be moved to the end
        self.assertTrue(id(g) == id(orig))

        # Adding another graph
        g += Graph.Full(3)
        self.assertTrue(g.vcount() == 11 and g.ecount() == 10
                        and g.clusters().membership == [0,0,0,0,0,0,1,2,3,3,3])

        # Adding two graphs
        g += [Graph.Full(3), Graph.Full(2)]
        self.assertTrue(g.vcount() == 16 and g.ecount() == 14
                        and g.clusters().membership == [0,0,0,0,0,0,1,2,3,3,3,4,4,4,5,5])

    def testAddition(self):
        g0 = Graph.Full(3)

        # Adding vertices
        g = g0+2
        self.assertTrue(g.vcount() == 5 and g.ecount() == 3
                        and g.clusters().membership == [0,0,0,1,2])
        g0 = g

        # Adding vertices by name
        g = g0+"spam"
        self.assertTrue(g.vcount() == 6 and g.ecount() == 3
                        and g.clusters().membership == [0,0,0,1,2,3])
        g0 = g

        # Adding a single edge
        g = g0+(2,3)
        self.assertTrue(g.vcount() == 6 and g.ecount() == 4
                        and g.clusters().membership == [0,0,0,0,1,2])
        g0 = g

        # Adding two edges
        g = g0+[(3, 4), (2, 4), (4, 5)]
        self.assertTrue(g.vcount() == 6 and g.ecount() == 7
                        and g.clusters().membership == [0]*6)
        g0 = g

        # Adding another graph
        g = g0+Graph.Full(3)
        self.assertTrue(g.vcount() == 9 and g.ecount() == 10
                        and g.clusters().membership == [0,0,0,0,0,0,1,1,1])

    def testInPlaceSubtraction(self):
        g = Graph.Full(8)
        orig = g

        # Deleting a vertex by vertex selector
        g -= 7
        self.assertTrue(g.vcount() == 7 and g.ecount() == 21
                        and g.clusters().membership == [0,0,0,0,0,0,0])

        # Deleting a vertex
        g -= g.vs[6]
        self.assertTrue(g.vcount() == 6 and g.ecount() == 15
                        and g.clusters().membership == [0,0,0,0,0,0])

        # Deleting two vertices
        g -= [4, 5]
        self.assertTrue(g.vcount() == 4 and g.ecount() == 6
                        and g.clusters().membership == [0,0,0,0])

        # Deleting an edge
        g -= (1, 2)
        self.assertTrue(g.vcount() == 4 and g.ecount() == 5
                        and g.clusters().membership == [0,0,0,0])

        # Deleting three more edges
        g -= [(1, 3), (0, 2), (0, 3)]
        self.assertTrue(g.vcount() == 4 and g.ecount() == 2
                        and g.clusters().membership == [0,0,1,1])
        
        # Did we really use the original graph so far?
        self.assertTrue(id(g) == id(orig))

        # Subtracting a graph
        g2 = Graph.Tree(3, 2)
        g -= g2
        self.assertTrue(g.vcount() == 4 and g.ecount() == 1
                        and g.clusters().membership == [0,1,2,2])

    def testNonzero(self):
        self.assertTrue(Graph(1))
        self.assertFalse(Graph(0))

    def testLength(self):
        self.assertRaises(TypeError, len, Graph(15))
        self.assertTrue(len(Graph(15).vs) == 15)
        self.assertTrue(len(Graph.Full(5).es) == 10)

    def testSimplify(self):
        el = [(0,1), (1,0), (1,2), (2,3), (2,3), (2,3), (3,3)] 
        g = Graph(el)
        g.es["weight"] = [1, 2, 3, 4, 5, 6, 7]

        g2 = g.copy()
        g2.simplify()
        self.assertTrue(g2.vcount() == g.vcount())
        self.assertTrue(g2.ecount() == 3)

        g2 = g.copy()
        g2.simplify(loops=False)
        self.assertTrue(g2.vcount() == g.vcount())
        self.assertTrue(g2.ecount() == 4)

        g2 = g.copy()
        g2.simplify(multiple=False)
        self.assertTrue(g2.vcount() == g.vcount())
        self.assertTrue(g2.ecount() == g.ecount() - 1)

    def testContractVertices(self):
        g = Graph.Full(4) + Graph.Full(4) + [(0, 5), (1, 4)]

        g2 = g.copy()
        g2.contract_vertices([0, 1, 2, 3, 1, 0, 4, 5])
        self.assertEqual(g2.vcount(), 6)
        self.assertEqual(g2.ecount(), g.ecount())
        self.assertEqual(sorted(g2.get_edgelist()),
                [(0, 0), (0, 1), (0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
                 (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (4, 5)])

        g2 = g.copy()
        g2.contract_vertices([0, 1, 2, 3, 1, 0, 6, 7])
        self.assertEqual(g2.vcount(), 8)
        self.assertEqual(g2.ecount(), g.ecount())
        self.assertEqual(sorted(g2.get_edgelist()),
                [(0, 0), (0, 1), (0, 1), (0, 2), (0, 3), (0, 6), (0, 7),
                 (1, 1), (1, 2), (1, 3), (1, 6), (1, 7), (2, 3), (6, 7)])

        g2 = Graph(10)
        g2.contract_vertices([0, 0, 1, 1, 2, 2, 3, 3, 4, 4])
        self.assertEqual(g2.vcount(), 5)
        self.assertEqual(g2.ecount(), 0)


def suite():
    operator_suite = unittest.makeSuite(OperatorTests)
    return unittest.TestSuite([operator_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

