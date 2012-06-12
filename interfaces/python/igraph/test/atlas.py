from __future__ import division

import warnings
import unittest
from igraph import *

class TestBase(unittest.TestCase):
    def testPageRank(self):
        for idx, g in enumerate(self.__class__.graphs):
            try:
                pr = g.pagerank()
            except Exception, ex:
                self.assertTrue(False, msg="PageRank calculation threw exception for graph #%d: %s" % (idx, ex))
                raise

            if g.vcount() == 0:
                self.assertEqual([], pr)
                continue

            self.assertAlmostEqual(1.0, sum(pr), places=5, \
                    msg="PageRank sum is not 1.0 for graph #%d (%r)" % (idx, pr))
            self.assertTrue(min(pr) >= 0, \
                    msg="Minimum PageRank is less than 0 for graph #%d (%r)" % (idx, pr))

    def testEigenvectorCentrality(self):
        # Temporarily turn off the warning handler because g.evcent() will print
        # a warning for DAGs
        warnings.simplefilter("ignore")

        try:
            for idx, g in enumerate(self.__class__.graphs):
                try:
                    ec, eval = g.evcent(return_eigenvalue=True)
                except Exception, ex:
                    self.assertTrue(False, msg="Eigenvector centrality threw exception for graph #%d: %s" % (idx, ex))
                    raise

                if g.vcount() == 0:
                    self.assertEqual([], ec)
                    continue

                if not g.is_connected():
                    # Skip disconnected graphs; this will be fixed in igraph 0.7
                    continue

                n = g.vcount()
                if abs(eval) < 1e-4:
                    self.assertTrue(min(ec) >= -1e-10,
                            msg="Minimum eigenvector centrality is smaller than 0 for graph #%d" % idx)
                    self.assertTrue(max(ec) <= 1,
                            msg="Maximum eigenvector centrality is greater than 1 for graph #%d" % idx)
                    continue

                self.assertAlmostEqual(max(ec), 1, places=7, \
                        msg="Maximum eigenvector centrality is %r (not 1) for graph #%d (%r)" % \
                        (max(ec), idx, ec))
                self.assertTrue(min(ec) >= 0, \
                        msg="Minimum eigenvector centrality is less than 0 for graph #%d" % idx)

                ec2 = [sum(ec[u.index] for u in v.predecessors()) for v in g.vs]
                for i in xrange(n):
                    self.assertAlmostEqual(ec[i] * eval, ec2[i], places=7, \
                            msg="Eigenvector centrality in graph #%d seems to be invalid "\
                            "for vertex %d" % (idx, i))
        finally:
            # Reset the warning handler
            warnings.resetwarnings()

    def testHubScore(self):
        for idx, g in enumerate(self.__class__.graphs):
            sc = g.hub_score()
            if g.vcount() == 0:
                self.assertEqual([], sc)
                continue

            self.assertAlmostEqual(max(sc), 1, places=7, \
                    msg="Maximum authority score is not 1 for graph #%d" % idx)
            self.assertTrue(min(sc) >= 0, \
                    msg="Minimum hub score is less than 0 for graph #%d" % idx)

    def testAuthorityScore(self):
        for idx, g in enumerate(self.__class__.graphs):
            sc = g.authority_score()
            if g.vcount() == 0:
                self.assertEqual([], sc)
                continue

            self.assertAlmostEqual(max(sc), 1, places=7, \
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

