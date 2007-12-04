import unittest
from igraph import *

class IsomorphismTests(unittest.TestCase):
    def testIsomorphic(self):
        g1 = Graph(8, [(0, 4), (0, 5), (0, 6), \
                       (1, 4), (1, 5), (1, 7), \
                       (2, 4), (2, 6), (2, 7), \
                       (3, 5), (3, 6), (3, 7)])
        g2 = Graph(8, [(0, 1), (0, 3), (0, 4), \
                       (2, 3), (2, 1), (2, 6), \
                       (5, 1), (5, 4), (5, 6), \
                       (7, 3), (7, 6), (7, 4)])
        self.failUnless(g1.isomorphic(g2))
        self.failUnless(g2.isomorphic(g1, return_mapping_21=True) == (True, None, [0, 2, 5, 7, 1, 3, 4, 6]))

    def testCountIsomorphisms(self):
        g = Graph.Full(4)
        self.failUnless(g.count_automorphisms() == 24)
        g = Graph(6, [(0,1), (2,3), (4,5), (0,2), (2,4), (1,3), (3,5)])
        self.failUnless(g.count_automorphisms() == 4)

    def testGetIsomorphisms(self):
        g = Graph(6, [(0,1), (2,3), (4,5), (0,2), (2,4), (1,3), (3,5)])
        maps = g.get_automorphisms()
        expected_maps = [[0,1,2,3,4,5], [1,0,3,2,5,4], [4,5,2,3,0,1], [5,4,3,2,1,0]]
        self.failUnless(maps == expected_maps)

    def testSubisomorphic(self):
        g = Graph.Lattice([3,3], circular=False)
        g2 = Graph([(0,1), (1,2), (1,3)])
        self.failUnless(g.subisomorphic(g2))
        self.failUnless(not g2.subisomorphic(g))

    def testCountSubisomorphisms(self):
        g = Graph.Lattice([3,3], circular=False)
        g2 = Graph.Lattice([2,2], circular=False)
        self.failUnless(g.count_subisomorphisms(g2) == 4*4*2)
        self.failUnless(g2.count_subisomorphisms(g) == 0)

def suite():
    isomorphism_suite = unittest.makeSuite(IsomorphismTests)
    return unittest.TestSuite([isomorphism_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

