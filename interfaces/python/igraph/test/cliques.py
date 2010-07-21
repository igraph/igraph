import unittest
from igraph import *

class CliqueTests(unittest.TestCase):
    def setUp(self):
        self.g=Graph.Full(6)
        self.g.delete_edges([(0, 1), (0, 2), (3, 5)])

    def testCliques(self):
        tests = {(4, -1): [[1, 2, 3, 4], [1, 2, 4, 5]],
                 (2, 2): [[0, 3], [0, 4], [0, 5],
                          [1, 2], [1, 3], [1, 4], [1, 5],
                          [2, 3], [2, 4], [2, 5], [3, 4], [4, 5]],
                 (-1, -1): [[0], [1], [2], [3], [4], [5],
                            [0, 3], [0, 4], [0, 5],
                            [1, 2], [1, 3], [1, 4], [1, 5],
                            [2, 3], [2, 4], [2, 5], [3, 4], [4, 5],
                            [0, 3, 4], [0, 4, 5],
                            [1, 2, 3], [1, 2, 4], [1, 2, 5],
                            [1, 3, 4], [1, 4, 5], [2, 3, 4], [2, 4, 5],
                            [1, 2, 3, 4], [1, 2, 4, 5]]}
        for (lo, hi), exp in tests.iteritems():
            self.assertEqual(sorted(exp), sorted(map(sorted, self.g.cliques(lo, hi))))

    def testLargestCliques(self):
        self.assertEqual(sorted(map(sorted, self.g.largest_cliques())),
                         [[1, 2, 3, 4], [1, 2, 4, 5]])

    def testMaximalCliques(self):
        self.assertEqual(sorted(map(sorted, self.g.maximal_cliques())),
                         [[0, 3, 4], [0, 4, 5],
                          [1, 2, 3, 4], [1, 2, 4, 5]])
        self.assertEqual(sorted(map(sorted, self.g.maximal_cliques(min=4))),
                         [[1, 2, 3, 4], [1, 2, 4, 5]])
        self.assertEqual(sorted(map(sorted, self.g.maximal_cliques(max=3))),
                         [[0, 3, 4], [0, 4, 5]])

    def testCliqueNumber(self):
        self.assertEqual(self.g.clique_number(), 4)
        self.assertEqual(self.g.omega(), 4)

class IndependentVertexSetTests(unittest.TestCase):
    def setUp(self):
        self.g1=Graph.Tree(5, 2, TREE_UNDIRECTED)
        self.g2=Graph.Tree(10, 2, TREE_UNDIRECTED)

    def testIndependentVertexSets(self):
        tests = {(4, -1): [],
                 (2, 2): [(0, 3), (0, 4), (1, 2), (2, 3), (2, 4), (3, 4)],
                 (-1, -1): [(0,), (1,), (2,), (3,), (4,),
                            (0, 3), (0, 4), (1, 2), (2, 3), (2, 4),
                            (3, 4), (0, 3, 4), (2, 3, 4)]}
        for (lo, hi), exp in tests.iteritems():
            self.assertEqual(exp, self.g1.independent_vertex_sets(lo, hi))

    def testLargestIndependentVertexSets(self):
        self.assertEqual(self.g1.largest_independent_vertex_sets(),
                         [(0, 3, 4), (2, 3, 4)])

    def testMaximalIndependentVertexSets(self):
        self.assertEqual(self.g2.maximal_independent_vertex_sets(),
                         [(0, 3, 4, 5, 6), (0, 3, 5, 6, 9),
                          (0, 4, 5, 6, 7, 8), (0, 5, 6, 7, 8, 9),
                          (1, 2, 7, 8, 9), (1, 5, 6, 7, 8, 9),
                          (2, 3, 4), (2, 3, 9), (2, 4, 7, 8)])

    def testIndependenceNumber(self):
        self.assertEqual(self.g2.independence_number(), 6)
        self.assertEqual(self.g2.alpha(), 6)


class MotifTests(unittest.TestCase):
    def setUp(self):
        self.g = Graph.Erdos_Renyi(100, 0.2, directed=True)

    def testDyads(self):
        """
        @note: this test is not exhaustive, it only checks whether the
          L{DyadCensus} objects "understand" attribute and item accessors
        """
        dc = self.g.dyad_census()
        accessors = ["mut", "mutual", "asym", "asymm", "asymmetric", "null"]
        for a in accessors:
            self.failUnless(isinstance(getattr(dc, a), int))
            self.failUnless(isinstance(dc[a], int))
        self.failUnless(isinstance(list(dc), list))
        self.failUnless(isinstance(tuple(dc), tuple))
        self.failUnless(len(list(dc)) == 3)
        self.failUnless(len(tuple(dc)) == 3)

    def testTriads(self):
        """
        @note: this test is not exhaustive, it only checks whether the
          L{TriadCensus} objects "understand" attribute and item accessors
        """
        tc = self.g.triad_census()
        accessors = ["003", "012", "021d", "030C"]
        for a in accessors:
            self.failUnless(isinstance(getattr(tc, "t"+a), int))
            self.failUnless(isinstance(tc[a], int))
        self.failUnless(isinstance(list(tc), list))
        self.failUnless(isinstance(tuple(tc), tuple))
        self.failUnless(len(list(tc)) == 16)
        self.failUnless(len(tuple(tc)) == 16)

class CliqueBenchmark(object):
    """This is a benchmark, not a real test case. You can run it
    using:

    >>> from igraph.test.clique import CliqueBenchmark
    >>> CliqueBenchmark().run()
    """

    def __init__(self):
        from time import time
        import gc
        self.time = time
        self.gc_collect = gc.collect

    def run(self):
        self.printIntro()
        self.testRandom()
        self.testMoonMoser()
        self.testGRG()

    def printIntro(self):
        print "n = number of vertices"
        print "#cliques = number of maximal cliques found"
        print "t1 = time required to determine the clique number"
        print "t2 = time required to determine and save all maximal cliques"
        print

    def timeit(self, g):
        start = self.time()
        omega = g.clique_number()
        mid = self.time()
        cl = g.maximal_cliques()
        end = self.time()
        self.gc_collect()
        return len(cl), mid-start, end-mid

    def testRandom(self):
        np = {100: [0.6, 0.7],
              300: [0.1, 0.2, 0.3, 0.4],
              500: [0.1, 0.2, 0.3],
              700: [0.1, 0.2],
              1000:[0.1, 0.2],
              10000: [0.001, 0.003, 0.005, 0.01, 0.02]}

        print
        print "Erdos-Renyi random graphs"
        print "       n        p #cliques        t1        t2"
        for n in sorted(np.keys()):
            for p in np[n]:
                g = Graph.Erdos_Renyi(n, p)
                result = self.timeit(g)
                print "%8d %8.3f %8d %8.4fs %8.4fs" % \
                    tuple([n, p] + list(result))

    def testMoonMoser(self):
        ns = [15, 27, 33]

        print
        print "Moon-Moser graphs"
        print "       n exp_clqs #cliques        t1        t2"
        for n in ns:
            n3 = n/3
            types = range(n3) * 3
            el = [(i, j) for i in range(n) for j in range(i+1,n) if types[i] != types[j]]
            g = Graph(n, el)
            result = self.timeit(g)
            print "%8d %8d %8d %8.4fs %8.4fs" % \
                tuple([n, (3**(n/3))] + list(result))

    def testGRG(self):
        ns = [100, 1000, 5000, 10000, 25000, 50000]

        print
        print "Geometric random graphs"
        print "       n        d #cliques        t1        t2"
        for n in ns:
            d = 2. / (n ** 0.5)
            g = Graph.GRG(n, d)
            result = self.timeit(g)
            print "%8d %8.3f %8d %8.4fs %8.4fs" % \
                tuple([n, d] + list(result))


def suite():
    clique_suite = unittest.makeSuite(CliqueTests)
    indvset_suite = unittest.makeSuite(IndependentVertexSetTests)
    motif_suite = unittest.makeSuite(MotifTests)
    return unittest.TestSuite([clique_suite, indvset_suite, motif_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

