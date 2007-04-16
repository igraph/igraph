import unittest
from igraph import *

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


class CommunityTests(unittest.TestCase):
    def testClauset(self):
        g = Graph.Full(5) + Graph.Full(5)
        g.add_edges([(0, 5)])
        cl = g.community_clauset()
        self.failUnless(cl.membership == [0,0,0,0,0,1,1,1,1,1])
        self.assertAlmostEqual(cl.q, 0.4523, places=3)
        self.failUnless(len(g.community_clauset(3)) == 3)
        
def suite():
    decomposition_suite = unittest.makeSuite(DecompositionTests)
    clustering_suite = unittest.makeSuite(ClusteringTests)
    community_suite = unittest.makeSuite(CommunityTests)
    return unittest.TestSuite([decomposition_suite, community_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

