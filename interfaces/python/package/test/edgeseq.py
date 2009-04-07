# vim:ts=4 sw=4 sts=4:
import unittest
from igraph import *

class EdgeSeqTests(unittest.TestCase):
    def setUp(self):
        self.g = Graph.Full(10)
        self.g.es["test"] = range(45)
    
    def testCreation(self):
        self.failUnless(len(EdgeSeq(self.g)) == 45)
        self.failUnless(len(EdgeSeq(self.g, 2)) == 1)
        self.failUnless(len(EdgeSeq(self.g, [1,2,3])) == 3)
        self.failUnless(EdgeSeq(self.g, [1,2,3]).indices == [1,2,3])
        self.assertRaises(ValueError, EdgeSeq, self.g, 112)
        self.assertRaises(ValueError, EdgeSeq, self.g, [112])
        self.failUnless(self.g.es.graph == self.g)

    def testPartialAttributeAssignment(self):
        only_even = self.g.es.select(lambda e: (e.index % 2 == 0))
        
        only_even["test"] = [0]*len(only_even)
        expected = [[0,i][i % 2] for i in xrange(self.g.ecount())]
        self.failUnless(self.g.es["test"] == expected)
        
        only_even["test2"] = range(23)
        expected = [[i/2, None][i % 2] for i in xrange(self.g.ecount())]
        self.failUnless(self.g.es["test2"] == expected)
        
        self.assertRaises(ValueError, only_even.__setitem__, "test2", range(6))

    def testAllSequence(self):
        self.failUnless(len(self.g.es) == 45)
        self.failUnless(self.g.es["test"] == range(45))

    def testEmptySequence(self):
        empty_es = self.g.es.select(None)
        self.failUnless(len(empty_es) == 0)
        self.assertRaises(IndexError, empty_es.__getitem__, 0)
        self.assertRaises(KeyError, empty_es.__getitem__, "nonexistent")
        self.failUnless(empty_es["test"] == [])
        empty_es = self.g.es[[]]
        self.failUnless(len(empty_es) == 0)
        empty_es = self.g.es[()]
        self.failUnless(len(empty_es) == 0)

    def testCallableFiltering(self):
        only_even = self.g.es.select(lambda e: (e.index % 2 == 0))
        self.failUnless(len(only_even) == 23)
        self.assertRaises(KeyError, only_even.__getitem__, "nonexistent")
        self.failUnless(only_even["test"] == [i*2 for i in xrange(23)])

    def testChainedCallableFiltering(self):
        only_div_six = self.g.es.select(lambda e: (e.index % 2 == 0),
          lambda e: (e.index % 3 == 0))
        self.failUnless(len(only_div_six) == 8)
        self.failUnless(only_div_six["test"] == [0, 6, 12, 18, 24, 30, 36, 42])

        only_div_six = self.g.es.select(lambda e: (e.index % 2 == 0)).select(\
          lambda e: (e.index % 3 == 0))
        self.failUnless(len(only_div_six) == 8)
        self.failUnless(only_div_six["test"] == [0, 6, 12, 18, 24, 30, 36, 42])

    def testIntegerFiltering(self):
        subset = self.g.es.select(2,3,4,2)
        self.failUnless(len(subset) == 4)
        self.failUnless(subset["test"] == [2,3,4,2])
        self.assertRaises(TypeError, self.g.es, "select", 2, 3, 4, 2, None)

        subset = self.g.es[2,3,4,2]
        self.failUnless(len(subset) == 4)
        self.failUnless(subset["test"] == [2,3,4,2])

    def testIterableFiltering(self):
        subset = self.g.es.select(xrange(5,8))
        self.failUnless(len(subset) == 3)
        self.failUnless(subset["test"] == [5,6,7])

    def testSliceFiltering(self):
        subset = self.g.es.select(slice(5, 8))
        self.failUnless(len(subset) == 3)
        self.failUnless(subset["test"] == [5,6,7])
        subset = self.g.es[40:56:2]
        self.failUnless(len(subset) == 3)
        self.failUnless(subset["test"] == [40,42,44])

    def testKeywordFiltering(self):
        g = Graph.Barabasi(1000, 2)
        g.es["betweenness"] = g.edge_betweenness()
        g.es["parity"] = [i % 2 for i in xrange(g.ecount())]
        self.failUnless(len(g.es(betweenness_gt=10)) < 2000)
        self.failUnless(len(g.es(betweenness_gt=10, parity=0)) < 2000)

    def testGraphMethodProxying(self):
        g = Graph.Barabasi(100)
        es = g.es(1,3,5,7,9)
        ebs = g.edge_betweenness()
        self.failUnless([ebs[i] for i in [1,3,5,7,9]] == es.edge_betweenness())


def suite():
    es_suite = unittest.makeSuite(EdgeSeqTests)
    return unittest.TestSuite([es_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

