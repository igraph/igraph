# vim:ts=4 sw=4 sts=4:
import unittest
from igraph import *

class EdgeTests(unittest.TestCase):
    def setUp(self):
        self.g = Graph.Full(10)

    def testRepr(self):
        output = repr(self.g.es[0])
        self.assertEquals(output, "igraph.Edge(%r, 0, {})" % self.g)

        self.g.es["weight"] = range(10, 0, -1)
        output = repr(self.g.es[3])
        self.assertEquals(output, "igraph.Edge(%r, 3, {'weight': 7})" % self.g)


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
        expected = [[i//2, None][i % 2] for i in xrange(self.g.ecount())]
        self.failUnless(self.g.es["test2"] == expected)

    def testSequenceReusing(self):
        if "test" in self.g.edge_attributes(): del self.g.es["test"]

        self.g.es["test"] = ["A", "B", "C"]
        self.failUnless(self.g.es["test"] == ["A", "B", "C"]*15)
        self.g.es["test"] = "ABC"
        self.failUnless(self.g.es["test"] == ["ABC"] * 45)

        only_even = self.g.es.select(lambda e: (e.index % 2 == 0))
        only_even["test"] = ["D", "E"]
        expected = ["D", "ABC", "E", "ABC"] * 12
        expected = expected[0:45]
        self.failUnless(self.g.es["test"] == expected)
        del self.g.es["test"]
        only_even["test"] = ["D", "E"]
        expected = ["D", None, "E", None] * 12
        expected = expected[0:45]
        self.failUnless(self.g.es["test"] == expected)

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

    def testSourceTargetFiltering(self):
        g = Graph.Barabasi(1000, 2)
        es1 = set(e.source for e in g.es.select(_target_in = [2,4]))
        es2 = set(v1 for v1, v2 in g.get_edgelist() if v2 in [2, 4])
        self.failUnless(es1 == es2)

    def testWithinFiltering(self):
        g = Graph.Lattice([10, 10])
        vs = [0, 1, 2, 10, 11, 12, 20, 21, 22] 

        es1 = g.es.select(_within = vs)
        es2 = g.es.select(_within = VertexSeq(g, vs))

        for es in [es1, es2]:
            self.failUnless(len(es) == 12)
            self.failUnless(all(e.source in vs and e.target in vs for e in es))

    def testBetweenFiltering(self):
        g = Graph.Lattice([10, 10])
        vs1, vs2 = [10, 11, 12], [20, 21, 22] 

        es1 = g.es.select(_between = (vs1, vs2))
        es2 = g.es.select(_between = (VertexSeq(g, vs1), VertexSeq(g, vs2)))

        for es in [es1, es2]:
            self.failUnless(len(es) == 3)
            self.failUnless(all((e.source in vs1 and e.target in vs2) or \
                                (e.target in vs1 and e.source in vs2) for e in es))

    def testGraphMethodProxying(self):
        g = Graph.Barabasi(100)
        es = g.es(1,3,5,7,9)
        ebs = g.edge_betweenness()
        self.failUnless([ebs[i] for i in [1,3,5,7,9]] == es.edge_betweenness())


def suite():
    edge_suite = unittest.makeSuite(EdgeTests)
    es_suite = unittest.makeSuite(EdgeSeqTests)
    return unittest.TestSuite([edge_suite, es_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

