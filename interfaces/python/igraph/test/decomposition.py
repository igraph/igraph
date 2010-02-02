import unittest
import math

from igraph import *
try:
    set, frozenset
except NameError:
    import sets
    set, frozenset = sets.Set, sets.ImmutableSet


class SubgraphTests(unittest.TestCase):
    def testSubgraph(self):
        g = Graph.Lattice([10, 10], circular=False, mutual=False)
        g.vs["id"] = range(g.vcount())

        vs = [0, 1, 2, 10, 11, 12, 20, 21, 22]
        sg = g.subgraph(vs)

        self.failUnless(sg.isomorphic(Graph.Lattice([3, 3], circular=False, mutual=False)))
        self.failUnless(sg.vs["id"] == vs)

    def testSubgraphEdges(self):
        g = Graph.Lattice([10, 10], circular=False, mutual=False)
        g.es["id"] = range(g.ecount())

        es = [0, 1, 2, 5, 20, 21, 22, 24, 38, 40]
        sg = g.subgraph_edges(es)
        exp = Graph.Lattice([3, 3], circular=False, mutual=False)
        exp.delete_edges([7, 8])

        self.failUnless(sg.isomorphic(exp))
        self.failUnless(sg.es["id"] == es)


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
        self.failUnless(self.cl.membership == [0,0,0,1,1,2,1,1,4,4])

    def testClusteringSizes(self):
        self.failUnless(self.cl.sizes() == [3, 4, 1, 0, 2])
        self.failUnless(self.cl.sizes(2, 4, 1) == [1, 2, 4])
        self.failUnless(self.cl.size(2) == 1)

    def testClusteringHistogram(self):
        self.failUnless(isinstance(self.cl.size_histogram(), Histogram))


class VertexClusteringTests(unittest.TestCase):
    def setUp(self):
        self.graph = Graph.Full(10)
        self.graph.vs["string"] = list("aaabbcccab")
        self.graph.vs["int"] = [17, 41, 23, 25, 64, 33, 3, 24, 47, 15]

    def testFromStringAttribute(self):
        cl = VertexClustering.FromAttribute(self.graph, "string")
        self.failUnless(cl.membership == [0,0,0,1,1,2,2,2,0,1])

    def testFromIntAttribute(self):
        cl = VertexClustering.FromAttribute(self.graph, "int")
        self.failUnless(cl.membership == list(range(10)))
        cl = VertexClustering.FromAttribute(self.graph, "int", 15)
        self.failUnless(cl.membership == [0, 1, 0, 0, 2, 1, 3, 0, 4, 0])
        cl = VertexClustering.FromAttribute(self.graph, "int", [10, 20, 30])
        self.failUnless(cl.membership == [0, 1, 2, 2, 1, 1, 3, 2, 1, 0])


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

    def testMultilevel(self):
        # Example graph from the paper
        g = Graph(16)
        g += [(0,2), (0,3), (0,4), (0,5),
              (1,2), (1,4), (1,7), (2,4), (2,5), (2,6),
              (3,7), (4,10), (5,7), (5,11), (6,7), (6,11),
              (8,9), (8,10), (8,11), (8,14), (8,15),
              (9,12), (9,14), (10,11), (10,12), (10,13),
              (10,14), (11,13)]
        cls = g.community_multilevel(return_levels=True)
        self.failUnless(len(cls) == 2)
        self.failUnless(cls[0].membership == [0,0,0,1,0,0,1,1,2,2,2,3,2,3,2,2])
        self.failUnless(cls[1].membership == [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1])
        self.assertAlmostEquals(cls[0].q, 0.346301, 5)
        self.assertAlmostEquals(cls[1].q, 0.392219, 5)

    
    def testWalktrap(self):
        g = Graph.Full(5) + Graph.Full(5) + Graph.Full(5)
        g += [(0,5), (5,10), (10, 0)]
        cl = g.community_walktrap()
        self.failUnless(cl.membership == [0,0,0,0,0,1,1,1,1,1,2,2,2,2,2])
        cl = g.community_walktrap(steps=3)
        self.failUnless(cl.membership == [0,0,0,0,0,1,1,1,1,1,2,2,2,2,2])
        
class ComparisonTests(unittest.TestCase):
    def testCompareVI(self):
        l1 = Clustering([1, 1, 1, 2, 2, 2])
        l2 = Clustering([1, 1, 2, 2, 3, 3])
        self.assertAlmostEqual(compare_communities(l1, l2), 0.8675, places=3)
        l1 = [1, 1, 1, 1, 1, 1]
        l2 = [1, 2, 3, 5, 6, 7]
        self.assertAlmostEqual(compare_communities(l1, l2), math.log(6), places=3)
        self.failUnless(compare_communities(l1, l1) == 0.)

def suite():
    decomposition_suite = unittest.makeSuite(DecompositionTests)
    clustering_suite = unittest.makeSuite(ClusteringTests)
    vertex_clustering_suite = unittest.makeSuite(VertexClusteringTests)
    overlapping_clustering_suite = unittest.makeSuite(OverlappingClusteringTests)
    community_suite = unittest.makeSuite(CommunityTests)
    comparison_suite = unittest.makeSuite(ComparisonTests)
    return unittest.TestSuite([decomposition_suite, clustering_suite, \
            vertex_clustering_suite, overlapping_clustering_suite, \
            community_suite, comparison_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

