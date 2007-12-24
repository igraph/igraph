import unittest
from igraph import *

class OperatorTests(unittest.TestCase):
    def testMultiplication(self):
        g = Graph.Full(3)*3
        self.failUnless(g.vcount() == 9 and g.ecount() == 9
                        and g.clusters().membership == [0,0,0,1,1,1,2,2,2])

    def testIntersection(self):
        g = Graph.Tree(7, 2) & Graph.Lattice([7])
        self.failUnless(g.get_edgelist() == [(0, 1)])

    def testUnion(self):
        g = Graph.Tree(7, 2) | Graph.Lattice([7])
        self.failUnless(g.vcount() == 7 and g.ecount() == 12)

    def testInPlaceAddition(self):
        g = Graph.Full(3)
        orig = g

        # Adding vertices
        g += 2
        self.failUnless(g.vcount() == 5 and g.ecount() == 3
                        and g.clusters().membership == [0,0,0,1,2])

        # Adding a single edge
        g += (2, 3)
        self.failUnless(g.vcount() == 5 and g.ecount() == 4
                        and g.clusters().membership == [0,0,0,0,1])

        # Adding two edges
        g += [(3, 4), (2, 4)]
        self.failUnless(g.vcount() == 5 and g.ecount() == 6
                        and g.clusters().membership == [0]*5)

        # Did we really use the original graph so far?
        # TODO: disjoint union should be modified so that this assertion
        # could be moved to the end
        self.assert_(id(g) == id(orig))

        # Adding another graph
        g += Graph.Full(3)
        self.failUnless(g.vcount() == 8 and g.ecount() == 9
                        and g.clusters().membership == [0,0,0,0,0,1,1,1])

        # Adding two graphs
        g += [Graph.Full(3), Graph.Full(2)]
        self.failUnless(g.vcount() == 13 and g.ecount() == 13
                        and g.clusters().membership == [0,0,0,0,0,1,1,1,2,2,2,3,3])

    def testAddition(self):
        g0 = Graph.Full(3)

        # Adding vertices
        g = g0+2
        self.failUnless(g.vcount() == 5 and g.ecount() == 3
                        and g.clusters().membership == [0,0,0,1,2])
        g0 = g

        # Adding a single edge
        g = g0+(2,3)
        self.failUnless(g.vcount() == 5 and g.ecount() == 4
                        and g.clusters().membership == [0,0,0,0,1])
        g0 = g

        # Adding two edges
        g = g0+[(3, 4), (2, 4)]
        self.failUnless(g.vcount() == 5 and g.ecount() == 6
                        and g.clusters().membership == [0]*5)
        g0 = g

        # Adding another graph
        g = g0+Graph.Full(3)
        self.failUnless(g.vcount() == 8 and g.ecount() == 9
                        and g.clusters().membership == [0,0,0,0,0,1,1,1])

    def testInPlaceSubtraction(self):
        g = Graph.Full(8)
        orig = g

        # Deleting a vertex by vertex selector
        g -= 7
        self.failUnless(g.vcount() == 7 and g.ecount() == 21
                        and g.clusters().membership == [0,0,0,0,0,0,0])

        # Deleting a vertex
        g -= g.vs[6]
        self.failUnless(g.vcount() == 6 and g.ecount() == 15
                        and g.clusters().membership == [0,0,0,0,0,0])

        # Deleting two vertices
        g -= [4, 5]
        self.failUnless(g.vcount() == 4 and g.ecount() == 6
                        and g.clusters().membership == [0,0,0,0])

        # Deleting an edge
        g -= (1, 2)
        self.failUnless(g.vcount() == 4 and g.ecount() == 5
                        and g.clusters().membership == [0,0,0,0])

        # Deleting three more edges
        g -= [(1, 3), (0, 2), (0, 3)]
        self.failUnless(g.vcount() == 4 and g.ecount() == 2
                        and g.clusters().membership == [0,0,1,1])
        
        # Did we really use the original graph so far?
        self.assert_(id(g) == id(orig))

        # Subtracting a graph
        g2 = Graph.Tree(3, 2)
        g -= g2
        self.failUnless(g.vcount() == 4 and g.ecount() == 1
                        and g.clusters().membership == [0,1,2,2])


def suite():
    operator_suite = unittest.makeSuite(OperatorTests)
    return unittest.TestSuite([operator_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

