import unittest
from igraph import *
try:
    set, frozenset
except NameError:
    import sets
    set, frozenset = sets.Set, sets.ImmutableSet

class DecompositionTests(unittest.TestCase):
    def testKCores(self):
        g = Graph(11, [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3),
                       (2,4), (2,5), (3,6), (3,7), (1,7), (7,8),
                       (1,9), (1,10), (9,10)])
        self.failUnless(g.coreness() == [3,3,3,3,1,1,1,2,1,2,2])
        self.failUnless(g.shell_index() == g.coreness())

        l=g.k_core(3).get_edgelist()
        l.sort()
        self.failUnless(l == [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)])

class ClusteringTests(unittest.TestCase):
    def setUp(self):
        self.cl = Clustering([0,0,0,1,1,2,1,1,4,4])

    def testClusteringIndex(self):
        self.failUnless(self.cl[0] == [0, 1, 2])
        self.failUnless(self.cl[1] == [3, 4, 6, 7])
        self.failUnless(self.cl[2] == [5])
        self.failUnless(self.cl[3] == [])
        self.failUnless(self.cl[4] == [8, 9])

    def testClusteringLength(self):
        self.failUnless(len(self.cl) == 5)

    def testClusteringMembership(self):
        self.failUnless(membership == [0,0,0,1,1,2,1,1,4,4])

    def testClusteringSizes(self):
        self.failUnless(self.cl.sizes() == [3, 4, 1, 0, 2])
        self.failUnless(self.cl.sizes(2, 4, 1) == [1, 2, 4])
        self.failUnless(self.cl.size(2) == 1)

    def testClusteringHistogram(self):
        self.failUnless(isinstance(self.cl.size_histogram(), Histogram))

class OverlappingClusteringTests(unittest.TestCase):
    def setUp(self):
        self.cl = OverlappingClustering([set([0]), set([0]), set([0]), set([0,1]), set([1]), set([1]), set([1]), set([]), set([3]), set([3,1])])

    def testClusteringIndex(self):
        self.failUnless(self.cl[0] == [0, 1, 2, 3])
        self.failUnless(self.cl[1] == [3, 4, 5, 6, 9])
        self.failUnless(self.cl[2] == [])
        self.failUnless(self.cl[3] == [8, 9])

    def testClusteringLength(self):
        self.failUnless(len(self.cl) == 4)

    def testClusteringSizes(self):
        self.failUnless(self.cl.sizes() == [4, 5, 0, 2])
        self.failUnless(self.cl.sizes(1, 3, 0) == [5, 2, 4])
        self.failUnless(self.cl.size(1) == 5)
        self.failUnless(self.cl.size(2) == 0)

    def testClusteringHistogram(self):
        self.failUnless(isinstance(self.cl.size_histogram(), Histogram))


class CommunityTests(unittest.TestCase):
    def testClauset(self):
        g = Graph.Full(5) + Graph.Full(5)
        g.add_edges([(0, 5)])
        cl = g.community_fastgreedy()
        self.failUnless(cl.membership == [0,0,0,0,0,1,1,1,1,1])
        self.assertAlmostEqual(cl.q, 0.4523, places=3)

        cl2 = OverlappingVertexClustering(g, [set([x]) for x in cl.membership])
        self.assertAlmostEqual(cl.q, cl2.q, places=3)

        g = Graph.Full(4) + Graph.Full(2)
        g.add_edges([(3,4)])
        weights = [1, 1, 1, 1, 1, 1, 10, 10]
        cl = g.community_fastgreedy(weights)
        self.failUnless(cl.membership == [0, 0, 0, 1, 1, 1])
        self.assertAlmostEqual(cl.q, 0.1708, places=3)
        
        g.es["weight"] = [3] * g.ecount()
        cl = g.community_fastgreedy("weight")
        self.failUnless(cl.membership == [0, 0, 0, 0, 1, 1])
        self.assertAlmostEqual(cl.q, 0.1796, places=3)

    def testEigenvector(self):
        g = Graph.Full(5) + Graph.Full(5)
        g.add_edges([(0, 5)])
        cl = g.community_leading_eigenvector_naive()
        self.failUnless(cl.membership == [0,0,0,0,0,1,1,1,1,1])
        self.assertAlmostEqual(cl.q, 0.4523, places=3)
        cl = g.community_leading_eigenvector(2)
        self.failUnless(cl.membership == [0,0,0,0,0,1,1,1,1,1])
        self.assertAlmostEqual(cl.q, 0.4523, places=3)

    def testLabelPropagation(self):
        # Nothing to test there really, since the algorithm
        # is pretty nondeterministic. We just do a quick smoke
        # test.
        g = Graph.GRG(100, 0.2)
        cl = g.community_label_propagation()
        g = Graph([(0,1),(1,2),(2,3)])
        g.es["weight"] = [2, 1, 2]
        g.vs["initial"] = [0, -1, -1, 1]
        cl = g.community_label_propagation("weight", "initial", [1,0,0,1])
        self.failUnless(cl.membership == [0, 0, 1, 1])
        cl = g.community_label_propagation(initial="initial", fixed=[1,0,0,1])
        self.failUnless(cl.membership == [0, 0, 1, 1] or \
                        cl.membership == [0, 1, 1, 1] or \
                        cl.membership == [0, 0, 0, 1])

    def testWalktrap(self):
        g = Graph.Full(5) + Graph.Full(5) + Graph.Full(5)
        g += [(0,5), (5,10), (10, 0)]
        cl = g.community_walktrap()
        self.failUnless(cl.membership == [0,0,0,0,0,1,1,1,1,1,2,2,2,2,2])
        cl = g.community_walktrap(steps=3)
        self.failUnless(cl.membership == [0,0,0,0,0,1,1,1,1,1,2,2,2,2,2])
        
def suite():
    decomposition_suite = unittest.makeSuite(DecompositionTests)
    clustering_suite = unittest.makeSuite(ClusteringTests)
    overlapping_clustering_suite = unittest.makeSuite(OverlappingClusteringTests)
    community_suite = unittest.makeSuite(CommunityTests)
    return unittest.TestSuite([decomposition_suite, community_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

