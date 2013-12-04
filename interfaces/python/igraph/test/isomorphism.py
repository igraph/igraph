import unittest
from igraph import *
from itertools import permutations
from random import shuffle

def node_compat(g1, g2, v1, v2):
    """Node compatibility function for isomorphism tests"""
    return g1.vs[v1]["color"] == g2.vs[v2]["color"]

def edge_compat(g1, g2, e1, e2):
    """Edge compatibility function for isomorphism tests"""
    return g1.es[e1]["color"] == g2.es[e2]["color"]

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
        self.assertTrue(g1.isomorphic(g2))
        self.assertTrue(g2.isomorphic_vf2(g1, return_mapping_21=True) \
          == (True, None, [0, 2, 5, 7, 1, 3, 4, 6]))
        self.assertTrue(g2.isomorphic_bliss(g1, return_mapping_21=True, sh2="fl")\
          == (True, None, [0, 2, 5, 7, 1, 3, 4, 6]))
        self.assertRaises(ValueError, g2.isomorphic_bliss, g1, sh2="nonexistent")

        # Test the automorphy of g1
        self.assertTrue(g1.isomorphic())
        self.assertTrue(g1.isomorphic_vf2(return_mapping_21=True) \
          == (True, None, [0, 1, 2, 3, 4, 5, 6, 7]))

        # Test VF2 with colors
        self.assertTrue(g1.isomorphic_vf2(g2,
            color1=[0,1,0,1,0,1,0,1],
            color2=[0,0,1,1,0,0,1,1]))
        g1.vs["color"] = [0,1,0,1,0,1,0,1]
        g2.vs["color"] = [0,0,1,1,0,1,1,0]
        self.assertTrue(not g1.isomorphic_vf2(g2, "color", "color"))

        # Test VF2 with vertex and edge colors
        self.assertTrue(g1.isomorphic_vf2(g2,
            color1=[0,1,0,1,0,1,0,1],
            color2=[0,0,1,1,0,0,1,1]))
        g1.es["color"] = range(12)
        g2.es["color"] = [0]*6 + [1]*6
        self.assertTrue(not g1.isomorphic_vf2(g2, "color", "color", "color", "color"))

        # Test VF2 with node compatibility function
        g2.vs["color"] = [0,0,1,1,0,0,1,1]
        self.assertTrue(g1.isomorphic_vf2(g2, node_compat_fn=node_compat))
        g2.vs["color"] = [0,0,1,1,0,1,1,0]
        self.assertTrue(not g1.isomorphic_vf2(g2, node_compat_fn=node_compat))

        # Test VF2 with node edge compatibility function
        g2.vs["color"] = [0,0,1,1,0,0,1,1]
        g1.es["color"] = range(12)
        g2.es["color"] = [0]*6 + [1]*6
        self.assertTrue(not g1.isomorphic_vf2(g2, node_compat_fn=node_compat,
            edge_compat_fn=edge_compat))

    def testIsomorphicCallback(self):
        maps = []
        def callback(g1, g2, map1, map2):
            maps.append(map1)
            return True

        # Test VF2 callback
        g = Graph(6, [(0,1), (2,3), (4,5), (0,2), (2,4), (1,3), (3,5)])
        g.isomorphic_vf2(g, callback=callback)
        expected_maps = [[0,1,2,3,4,5], [1,0,3,2,5,4], [4,5,2,3,0,1], [5,4,3,2,1,0]]
        self.assertTrue(sorted(maps) == expected_maps)

        maps[:] = []
        g3 = Graph.Full(4)
        g3.vs["color"] = [0,1,1,0]
        g3.isomorphic_vf2(callback=callback, color1="color", color2="color")
        expected_maps = [[0,1,2,3], [0,2,1,3], [3,1,2,0], [3,2,1,0]]
        self.assertTrue(sorted(maps) == expected_maps)

    def testCountIsomorphisms(self):
        g = Graph.Full(4)
        self.assertTrue(g.count_automorphisms_vf2() == 24)
        g = Graph(6, [(0,1), (2,3), (4,5), (0,2), (2,4), (1,3), (3,5)])
        self.assertTrue(g.count_automorphisms_vf2() == 4)

        # Some more tests with colors
        g3 = Graph.Full(4)
        g3.vs["color"] = [0,1,1,0]
        self.assertTrue(g3.count_isomorphisms_vf2() == 24)
        self.assertTrue(g3.count_isomorphisms_vf2(color1="color", color2="color") == 4)
        self.assertTrue(g3.count_isomorphisms_vf2(color1=[0,1,2,0], color2=(0,1,2,0)) == 2)
        self.assertTrue(g3.count_isomorphisms_vf2(edge_color1=[0,1,0,0,0,1],
            edge_color2=[0,1,0,0,0,1]) == 2)

        # Test VF2 with node/edge compatibility function
        g3.vs["color"] = [0,1,1,0]
        self.assertTrue(g3.count_isomorphisms_vf2(node_compat_fn=node_compat) == 4)
        g3.vs["color"] = [0,1,2,0]
        self.assertTrue(g3.count_isomorphisms_vf2(node_compat_fn=node_compat) == 2)
        g3.es["color"] = [0,1,0,0,0,1]
        self.assertTrue(g3.count_isomorphisms_vf2(edge_compat_fn=edge_compat) == 2)

    def testGetIsomorphisms(self):
        g = Graph(6, [(0,1), (2,3), (4,5), (0,2), (2,4), (1,3), (3,5)])
        maps = g.get_automorphisms_vf2()
        expected_maps = [[0,1,2,3,4,5], [1,0,3,2,5,4], [4,5,2,3,0,1], [5,4,3,2,1,0]]
        self.assertTrue(maps == expected_maps)

        g3 = Graph.Full(4)
        g3.vs["color"] = [0,1,1,0]
        expected_maps = [[0,1,2,3], [0,2,1,3], [3,1,2,0], [3,2,1,0]]
        self.assertTrue(sorted(g3.get_automorphisms_vf2(color="color")) == expected_maps)

class SubisomorphismTests(unittest.TestCase):
    def testSubisomorphicLAD(self):
        g = Graph.Lattice([3,3], circular=False)
        g2 = Graph([(0,1), (1,2), (1,3)])
        g3 = g + [(0,4), (2,4), (6,4), (8,4), (3,1), (1,5), (5,7), (7,3)]

        self.assertTrue(g.subisomorphic_lad(g2))
        self.assertFalse(g2.subisomorphic_lad(g))

        # Test 'induced'
        self.assertFalse(g3.subisomorphic_lad(g, induced=True))
        self.assertTrue(g3.subisomorphic_lad(g, induced=False))
        self.assertTrue(g3.subisomorphic_lad(g))
        self.assertTrue(g3.subisomorphic_lad(g2, induced=True))
        self.assertTrue(g3.subisomorphic_lad(g2, induced=False))
        self.assertTrue(g3.subisomorphic_lad(g2))

        # Test with limited vertex matching
        domains = [[4], [0,1,2,3,5,6,7,8], [0,1,2,3,5,6,7,8], [0,1,2,3,5,6,7,8]]
        self.assertTrue(g.subisomorphic_lad(g2, domains=domains))
        domains = [[], [0,1,2,3,5,6,7,8], [0,1,2,3,5,6,7,8], [0,1,2,3,5,6,7,8]]
        self.assertTrue(not g.subisomorphic_lad(g2, domains=domains))

    def testGetSubisomorphismsLAD(self):
        g = Graph.Lattice([3,3], circular=False)
        g2 = Graph([(0,1), (1,2), (2,3), (3,0)])
        g3 = g + [(0,4), (2,4), (6,4), (8,4), (3,1), (1,5), (5,7), (7,3)]

        all_subiso = "0143 0341 1034 1254 1430 1452 2145 2541 3014 3410 3476 \
        3674 4103 4125 4301 4367 4521 4587 4763 4785 5214 5412 5478 5874 6347 \
        6743 7436 7458 7634 7854 8547 8745"
        all_subiso = sorted([int(x) for x in item] for item in all_subiso.split())

        self.assertEqual(all_subiso, sorted(g.get_subisomorphisms_lad(g2)))
        self.assertEqual([], sorted(g2.get_subisomorphisms_lad(g)))

        # Test 'induced'
        induced_subiso = "1375 1573 3751 5731 7513 7315 5137 3157"
        induced_subiso = sorted([int(x) for x in item] for item in induced_subiso.split())
        all_subiso_extra = sorted(all_subiso + induced_subiso)
        self.assertEqual(induced_subiso,
                sorted(g3.get_subisomorphisms_lad(g2, induced=True)))
        self.assertEqual([], g3.get_subisomorphisms_lad(g, induced=True))
        
        # Test with limited vertex matching
        limited_subiso = [iso for iso in all_subiso if iso[0] == 4]
        domains = [[4], [0,1,2,3,5,6,7,8], [0,1,2,3,5,6,7,8], [0,1,2,3,5,6,7,8]]
        self.assertEqual(limited_subiso,
                sorted(g.get_subisomorphisms_lad(g2, domains=domains)))
        domains = [[], [0,1,2,3,5,6,7,8], [0,1,2,3,5,6,7,8], [0,1,2,3,5,6,7,8]]
        self.assertEqual([], sorted(g.get_subisomorphisms_lad(g2, domains=domains)))

    def testSubisomorphicVF2(self):
        g = Graph.Lattice([3,3], circular=False)
        g2 = Graph([(0,1), (1,2), (1,3)])
        self.assertTrue(g.subisomorphic_vf2(g2))
        self.assertTrue(not g2.subisomorphic_vf2(g))

        # Test with vertex colors
        g.vs["color"] = [0,0,0,0,1,0,0,0,0]
        g2.vs["color"] = [1,0,0,0]
        self.assertTrue(g.subisomorphic_vf2(g2, node_compat_fn=node_compat))
        g2.vs["color"] = [2,0,0,0]
        self.assertTrue(not g.subisomorphic_vf2(g2, node_compat_fn=node_compat))

        # Test with edge colors
        g.es["color"] = [1] + [0]*(g.ecount()-1)
        g2.es["color"] = [1] + [0]*(g2.ecount()-1)
        self.assertTrue(g.subisomorphic_vf2(g2, edge_compat_fn=edge_compat))
        g2.es[0]["color"] = [2]
        self.assertTrue(not g.subisomorphic_vf2(g2, node_compat_fn=node_compat))

    def testCountSubisomorphisms(self):
        g = Graph.Lattice([3,3], circular=False)
        g2 = Graph.Lattice([2,2], circular=False)
        self.assertTrue(g.count_subisomorphisms_vf2(g2) == 4*4*2)
        self.assertTrue(g2.count_subisomorphisms_vf2(g) == 0)

        # Test with vertex colors
        g.vs["color"] = [0,0,0,0,1,0,0,0,0]
        g2.vs["color"] = [1,0,0,0]
        self.assertTrue(g.count_subisomorphisms_vf2(g2, "color", "color") == 4*2)
        self.assertTrue(g.count_subisomorphisms_vf2(g2, node_compat_fn=node_compat) == 4*2)

        # Test with edge colors
        g.es["color"] = [1] + [0]*(g.ecount()-1)
        g2.es["color"] = [1] + [0]*(g2.ecount()-1)
        self.assertTrue(g.count_subisomorphisms_vf2(g2, edge_color1="color", edge_color2="color") == 2)
        self.assertTrue(g.count_subisomorphisms_vf2(g2, edge_compat_fn=edge_compat) == 2)

class PermutationTests(unittest.TestCase):
    def testCanonicalPermutation(self):
        # Simple case: two ring graphs
        g1 = Graph(4, [(0, 1), (1, 2), (2, 3), (3, 0)])
        g2 = Graph(4, [(0, 1), (1, 3), (3, 2), (2, 0)])

        cp = g1.canonical_permutation()
        g3 = g1.permute_vertices(cp)

        cp = g2.canonical_permutation()
        g4 = g2.permute_vertices(cp)

        self.assertTrue(g3.vcount() == g4.vcount())
        self.assertTrue(sorted(g3.get_edgelist()) == sorted(g4.get_edgelist()))

        # More complicated one: small GRG, random permutation
        g = Graph.GRG(10, 0.5)
        perm = range(10)
        shuffle(perm)
        g2 = g.permute_vertices(perm)
        g3 = g.permute_vertices(g.canonical_permutation())
        g4 = g2.permute_vertices(g2.canonical_permutation())

        self.assertTrue(g3.vcount() == g4.vcount())
        self.assertTrue(sorted(g3.get_edgelist()) == sorted(g4.get_edgelist()))

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
        self.assertTrue(g3.vcount() == g2.vcount() and g3.ecount() == g2.ecount())
        self.assertTrue(set(g3.get_edgelist()) == set(g1.get_edgelist()))

def suite():
    isomorphism_suite = unittest.makeSuite(IsomorphismTests)
    subisomorphism_suite = unittest.makeSuite(SubisomorphismTests)
    permutation_suite = unittest.makeSuite(PermutationTests)
    return unittest.TestSuite([isomorphism_suite, subisomorphism_suite, \
            permutation_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

