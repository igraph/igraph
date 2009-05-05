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

        # Test the isomorphy of g1 and g2
        self.failUnless(g1.isomorphic(g2))
        self.failUnless(g2.isomorphic_vf2(g1, return_mapping_21=True) \
          == (True, None, [0, 2, 5, 7, 1, 3, 4, 6]))
        self.failUnless(g2.isomorphic_bliss(g1, return_mapping_21=True, sh2="fl")\
          == (True, None, [0, 2, 5, 7, 1, 3, 4, 6]))
        self.assertRaises(ValueError, g2.isomorphic_bliss, g1, sh2="nonexistent")

        # Test the automorphy of g1
        self.failUnless(g1.isomorphic())
        self.failUnless(g1.isomorphic_vf2(return_mapping_21=True) \
          == (True, None, [0, 1, 2, 3, 4, 5, 6, 7]))

        # Test VF2 with colors
        self.failUnless(g1.isomorphic_vf2(g2,
            color1=[0,1,0,1,0,1,0,1],
            color2=[0,0,1,1,0,0,1,1]))
        g1.vs["color"] = [0,1,0,1,0,1,0,1]
        g2.vs["color"] = [0,0,1,1,0,1,1,0]
        self.failUnless(not g1.isomorphic_vf2(g2, "color", "color"))

        # Test VF2 with vertex and edge colors
        self.failUnless(g1.isomorphic_vf2(g2,
            color1=[0,1,0,1,0,1,0,1],
            color2=[0,0,1,1,0,0,1,1]))
        g1.es["color"] = range(12)
        g2.es["color"] = [0]*6 + [1]*6
        self.failUnless(not g1.isomorphic_vf2(g2, "color", "color", "color", "color"))

    def testCountIsomorphisms(self):
        g = Graph.Full(4)
        self.failUnless(g.count_automorphisms_vf2() == 24)
        g = Graph(6, [(0,1), (2,3), (4,5), (0,2), (2,4), (1,3), (3,5)])
        self.failUnless(g.count_automorphisms_vf2() == 4)

        # Some more tests with colors
        g3 = Graph.Full(4)
        g3.vs["color"] = [0,1,1,0]
        self.failUnless(g3.count_isomorphisms_vf2() == 24)
        self.failUnless(g3.count_isomorphisms_vf2(color1="color", color2="color") == 4)
        self.failUnless(g3.count_isomorphisms_vf2(color1=[0,1,2,0], color2=(0,1,2,0)) == 2)
        self.failUnless(g3.count_isomorphisms_vf2(edge_color1=[0,1,0,0,0,1],
            edge_color2=[0,1,0,0,0,1]) == 2)

    def testGetIsomorphisms(self):
        g = Graph(6, [(0,1), (2,3), (4,5), (0,2), (2,4), (1,3), (3,5)])
        maps = g.get_automorphisms_vf2()
        expected_maps = [[0,1,2,3,4,5], [1,0,3,2,5,4], [4,5,2,3,0,1], [5,4,3,2,1,0]]
        self.failUnless(maps == expected_maps)

        g3 = Graph.Full(4)
        g3.vs["color"] = [0,1,1,0]
        expected_maps = [[0,1,2,3], [0,2,1,3], [3,1,2,0], [3,2,1,0]]
        self.failUnless(sorted(g3.get_automorphisms_vf2()), expected_maps)


    def testSubisomorphic(self):
        g = Graph.Lattice([3,3], circular=False)
        g2 = Graph([(0,1), (1,2), (1,3)])
        self.failUnless(g.subisomorphic_vf2(g2))
        self.failUnless(not g2.subisomorphic_vf2(g))

    def testCountSubisomorphisms(self):
        g = Graph.Lattice([3,3], circular=False)
        g2 = Graph.Lattice([2,2], circular=False)
        self.failUnless(g.count_subisomorphisms_vf2(g2) == 4*4*2)
        self.failUnless(g2.count_subisomorphisms_vf2(g) == 0)

        # Test with vertex colors
        g.vs["color"] = [0,0,0,0,1,0,0,0,0]
        g2.vs["color"] = [1,0,0,0]
        self.failUnless(g.count_subisomorphisms_vf2(g2, "color", "color") == 4*2)

        # Test with edge colors
        g.es["color"] = [1] + [0]*(g.ecount()-1)
        g2.es["color"] = [1] + [0]*(g2.ecount()-1)
        self.failUnless(g.count_subisomorphisms_vf2(g2, edge_color1="color", edge_color2="color") == 2)

    def testPermuteVertices(self):
        g1 = Graph(8, [(0, 4), (0, 5), (0, 6), \
                       (1, 4), (1, 5), (1, 7), \
                       (2, 4), (2, 6), (2, 7), \
                       (3, 5), (3, 6), (3, 7)])
        g2 = Graph(8, [(0, 1), (0, 3), (0, 4), \
                       (2, 3), (2, 1), (2, 6), \
                       (5, 1), (5, 4), (5, 6), \
                       (7, 3), (7, 6), (7, 4)])
        _, _, mapping = g1.isomorphic_vf2(g2, return_mapping_21=True)
        g3 = g2.permute_vertices(mapping)
        self.failUnless(g3.vcount() == g2.vcount() and g3.ecount() == g2.ecount())
        self.failUnless(set(g3.get_edgelist()) == set(g1.get_edgelist()))

def suite():
    isomorphism_suite = unittest.makeSuite(IsomorphismTests)
    return unittest.TestSuite([isomorphism_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

