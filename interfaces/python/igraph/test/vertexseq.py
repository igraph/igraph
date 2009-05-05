# vim:ts=4 sw=4 sts=4:
import unittest
from igraph import *

class VertexSeqTests(unittest.TestCase):
    def setUp(self):
        self.g = Graph.Full(10)
        self.g.vs["test"] = range(10)
    
    def testCreation(self):
        self.failUnless(len(VertexSeq(self.g)) == 10)
        self.failUnless(len(VertexSeq(self.g, 2)) == 1)
        self.failUnless(len(VertexSeq(self.g, [1,2,3])) == 3)
        self.failUnless(VertexSeq(self.g, [1,2,3]).indices == [1,2,3])
        self.assertRaises(ValueError, VertexSeq, self.g, 12)
        self.assertRaises(ValueError, VertexSeq, self.g, [12])
        self.failUnless(self.g.vs.graph == self.g)

    def testPartialAttributeAssignment(self):
        only_even = self.g.vs.select(lambda v: (v.index % 2 == 0))
        only_even["test"] = [0]*len(only_even)
        self.failUnless(self.g.vs["test"] == [0,1,0,3,0,5,0,7,0,9])
        only_even["test2"] = range(5)
        self.failUnless(self.g.vs["test2"] == [0,None,1,None,2,None,3,None,4,None])
        self.assertRaises(ValueError, only_even.__setitem__, "test2", range(6))

    def testAllSequence(self):
        self.failUnless(len(self.g.vs) == 10)
        self.failUnless(self.g.vs["test"] == range(10))

    def testEmptySequence(self):
        empty_vs = self.g.vs.select(None)
        self.failUnless(len(empty_vs) == 0)
        self.assertRaises(IndexError, empty_vs.__getitem__, 0)
        self.assertRaises(KeyError, empty_vs.__getitem__, "nonexistent")
        self.failUnless(empty_vs["test"] == [])
        empty_vs = self.g.vs[[]]
        self.failUnless(len(empty_vs) == 0)
        empty_vs = self.g.vs[()]
        self.failUnless(len(empty_vs) == 0)

    def testCallableFiltering(self):
        only_even = self.g.vs.select(lambda v: (v.index % 2 == 0))
        self.failUnless(len(only_even) == 5)
        self.assertRaises(KeyError, only_even.__getitem__, "nonexistent")
        self.failUnless(only_even["test"] == [0, 2, 4, 6, 8])

    def testChainedCallableFiltering(self):
        only_div_six = self.g.vs.select(lambda v: (v.index % 2 == 0),
          lambda v: (v.index % 3 == 0))
        self.failUnless(len(only_div_six) == 2)
        self.failUnless(only_div_six["test"] == [0, 6])

        only_div_six = self.g.vs.select(lambda v: (v.index % 2 == 0)).select(\
          lambda v: (v.index % 3 == 0))
        self.failUnless(len(only_div_six) == 2)
        self.failUnless(only_div_six["test"] == [0, 6])

    def testIntegerFiltering(self):
        subset = self.g.vs.select(2,3,4,2)
        self.failUnless(len(subset) == 4)
        self.failUnless(subset["test"] == [2,3,4,2])
        self.assertRaises(TypeError, self.g.vs, "select", 2, 3, 4, 2, None)

        subset = self.g.vs[2,3,4,2]
        self.failUnless(len(subset) == 4)
        self.failUnless(subset["test"] == [2,3,4,2])

    def testIterableFiltering(self):
        subset = self.g.vs.select(xrange(5,8))
        self.failUnless(len(subset) == 3)
        self.failUnless(subset["test"] == [5,6,7])

    def testSliceFiltering(self):
        subset = self.g.vs.select(slice(5, 8))
        self.failUnless(len(subset) == 3)
        self.failUnless(subset["test"] == [5,6,7])
        subset = self.g.vs[5:16:2]
        self.failUnless(len(subset) == 3)
        self.failUnless(subset["test"] == [5,7,9])

    def testKeywordFiltering(self):
        g = Graph.Barabasi(10000)
        g.vs["degree"] = g.degree()
        g.vs["parity"] = [i % 2 for i in xrange(g.vcount())]
        l = len(g.vs(degree_gt=30))
        self.failUnless(l < 1000)
        self.failUnless(len(g.vs(degree_gt=30, parity=0)) <= 500)
        del g.vs["degree"]
        self.failUnless(len(g.vs(_degree_gt=30)) == l)

    def testGraphMethodProxying(self):
        g = Graph.Barabasi(100)
        vs = g.vs(1,3,5,7,9)
        self.failUnless(vs.degree() == g.degree(vs))
        self.failUnless(g.degree(vs) == g.degree(vs.indices))

def suite():
    vs_suite = unittest.makeSuite(VertexSeqTests)
    return unittest.TestSuite([vs_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

