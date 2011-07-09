from __future__ import division

import unittest
from igraph import *

class TestBase(unittest.TestCase):
    def testPageRank(self):
        for idx, g in enumerate(self.__class__.graphs):
            pr = g.pagerank()
            if g.vcount() == 0:
                self.assertEquals([], pr)
                continue

            self.assertAlmostEquals(1.0, sum(pr), places=5, \
                    msg="PageRank sum is not 1.0 for graph #%d" % idx)
            self.assertTrue(min(pr) >= 0, \
                    msg="Minimum PageRank is less than 0 for graph #%d" % idx)

    def testEigenvectorCentrality(self):
        for idx, g in enumerate(self.__class__.graphs):
            try:
                ec, eval = g.evcent(return_eigenvalue=True)
            except Exception as ex:
                self.assertTrue(False, msg="Eigenvector centrality threw exception for graph #%d: %s" % (idx, ex))
                raise

            if g.vcount() == 0:
                self.assertEquals([], ec)
                continue

            n = g.vcount()
            if abs(eval) < 1e-4:
                self.assertTrue(min(ec) >= -1e-10,
                        msg="Minimum eigenvector centrality is smaller than 0 for graph #%d" % idx)
                self.assertTrue(max(ec) <= 1,
                        msg="Maximum eigenvector centrality is greater than 1 for graph #%d" % idx)
                continue

            self.assertAlmostEquals(max(ec), 1, places=7, \
                    msg="Maximum eigenvector centrality is not 1 for graph #%d" % idx)
            self.assertTrue(min(ec) >= 0, \
                    msg="Minimum eigenvector centrality is less than 0 for graph #%d" % idx)

            ec2 = [sum(ec[u.index] for u in v.predecessors()) for v in g.vs]
            for i in xrange(n):
                self.assertAlmostEquals(ec[i] * eval, ec2[i], places=7, \
                        msg="Eigenvector centrality in graph #%d seems to be invalid "\
                        "for vertex %d" % (idx, i))

    def testHubScore(self):
        for idx, g in enumerate(self.__class__.graphs):
            sc = g.hub_score()
            if g.vcount() == 0:
                self.assertEquals([], sc)
                continue

            self.assertAlmostEquals(max(sc), 1, places=7, \
                    msg="Maximum authority score is not 1 for graph #%d" % idx)
            self.assertTrue(min(sc) >= 0, \
                    msg="Minimum hub score is less than 0 for graph #%d" % idx)

    def testAuthorityScore(self):
        for idx, g in enumerate(self.__class__.graphs):
            sc = g.authority_score()
            if g.vcount() == 0:
                self.assertEquals([], sc)
                continue

            self.assertAlmostEquals(max(sc), 1, places=7, \
                    msg="Maximum authority score is not 1 for graph #%d" % idx)
            self.assertTrue(min(sc) >= 0, \
                    msg="Minimum authority score is less than 0 for graph #%d" % idx)

class GraphAtlasTests(TestBase):
    graphs = [Graph.Atlas(i) for i in xrange(1253)]

class IsoclassTests(TestBase):
    graphs = [Graph.Isoclass(3, i, directed=True) for i in xrange(16)] + \
             [Graph.Isoclass(4, i, directed=True) for i in xrange(218)]

def suite():
    atlas_suite = unittest.makeSuite(GraphAtlasTests)
    isoclass_suite = unittest.makeSuite(IsoclassTests)
    return unittest.TestSuite([atlas_suite, isoclass_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

