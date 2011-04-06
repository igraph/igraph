# vim:ts=4 sw=4 sts=4:
import unittest
from igraph import *

class VertexTests(unittest.TestCase):
    def setUp(self):
        self.g = Graph.Full(10)

    def testUpdateAttributes(self):
        v = self.g.vs[0]

        v.update_attributes(a=2)
        self.assertEquals(v["a"], 2)

        v.update_attributes([("a", 3), ("b", 4)], c=5, d=6)
        self.assertEquals(v.attributes(), dict(a=3, b=4, c=5, d=6))

        v.update_attributes(dict(b=44, c=55))
        self.assertEquals(v.attributes(), dict(a=3, b=44, c=55, d=6))


class VertexSeqTests(unittest.TestCase):
    def setUp(self):
        self.g = Graph.Full(10)
        self.g.vs["test"] = range(10)
        self.g.vs["name"] = list("ABCDEFGHIJ")
    
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

    def testSequenceReusing(self):
        if "test" in self.g.vertex_attributes():
            del self.g.vs["test"]

        self.g.vs["test"] = ["A", "B", "C"]
        self.failUnless(self.g.vs["test"] == ["A", "B", "C", "A", "B", "C", "A", "B", "C", "A"])
        self.g.vs["test"] = "ABC"
        self.failUnless(self.g.vs["test"] == ["ABC"] * 10)

        only_even = self.g.vs.select(lambda v: (v.index % 2 == 0))
        only_even["test"] = ["D", "E"]
        self.failUnless(self.g.vs["test"] == ["D", "ABC", "E", "ABC", "D", "ABC", "E", "ABC", "D", "ABC"])
        del self.g.vs["test"]
        only_even["test"] = ["D", "E"]
        self.failUnless(self.g.vs["test"] == ["D", None, "E", None, "D", None, "E", None, "D", None])



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

    def testCallableFilteringFind(self):
        vertex = self.g.vs.find(lambda v: (v.index % 2 == 1))
        self.failUnless(vertex.index == 1)
        self.assertRaises(IndexError, self.g.vs.find, lambda v: (v.index % 2 == 3))

    def testCallableFilteringSelect(self):
        only_even = self.g.vs.select(lambda v: (v.index % 2 == 0))
        self.failUnless(len(only_even) == 5)
        self.assertRaises(KeyError, only_even.__getitem__, "nonexistent")
        self.failUnless(only_even["test"] == [0, 2, 4, 6, 8])

    def testChainedCallableFilteringSelect(self):
        only_div_six = self.g.vs.select(lambda v: (v.index % 2 == 0),
          lambda v: (v.index % 3 == 0))
        self.failUnless(len(only_div_six) == 2)
        self.failUnless(only_div_six["test"] == [0, 6])

        only_div_six = self.g.vs.select(lambda v: (v.index % 2 == 0)).select(\
          lambda v: (v.index % 3 == 0))
        self.failUnless(len(only_div_six) == 2)
        self.failUnless(only_div_six["test"] == [0, 6])

    def testIntegerFilteringFind(self):
        self.assertEquals(self.g.vs.find(3).index, 3)
        self.assertEquals(self.g.vs.select(2,3,4,2).find(3).index, 2)
        self.assertRaises(IndexError, self.g.vs.find, 17)

    def testIntegerFilteringSelect(self):
        subset = self.g.vs.select(2,3,4,2)
        self.assertEquals(len(subset), 4)
        self.assertEquals(subset["test"], [2,3,4,2])
        self.assertRaises(TypeError, self.g.vs.select, 2, 3, 4, 2, None)

        subset = self.g.vs[2,3,4,2]
        self.failUnless(len(subset) == 4)
        self.failUnless(subset["test"] == [2,3,4,2])

    def testStringFilteringFind(self):
        self.assertEquals(self.g.vs.find("D").index, 3)
        self.assertEquals(self.g.vs.select(2,3,4,2).find("C").index, 2)
        self.assertRaises(ValueError, self.g.vs.select(2,3,4,2).find, "F")
        self.assertRaises(ValueError, self.g.vs.find, "NoSuchName")

    def testIterableFilteringSelect(self):
        subset = self.g.vs.select(xrange(5,8))
        self.failUnless(len(subset) == 3)
        self.failUnless(subset["test"] == [5,6,7])

    def testSliceFilteringSelect(self):
        subset = self.g.vs.select(slice(5, 8))
        self.failUnless(len(subset) == 3)
        self.failUnless(subset["test"] == [5,6,7])
        subset = self.g.vs[5:16:2]
        self.failUnless(len(subset) == 3)
        self.failUnless(subset["test"] == [5,7,9])

    def testKeywordFilteringSelect(self):
        g = Graph.Barabasi(10000)
        g.vs["degree"] = g.degree()
        g.vs["parity"] = [i % 2 for i in xrange(g.vcount())]
        l = len(g.vs(degree_gt=30))
        self.failUnless(l < 1000)
        self.failUnless(len(g.vs(degree_gt=30, parity=0)) <= 500)
        del g.vs["degree"]
        self.failUnless(len(g.vs(_degree_gt=30)) == l)

    def testIndexOutOfBoundsSelect(self):
        g = Graph.Full(3)
        self.assertRaises(ValueError, g.vs.select, 4)
        self.assertRaises(ValueError, g.vs.select, 4, 5)
        self.assertRaises(ValueError, g.vs.select, (4, 5))
        self.assertRaises(ValueError, g.vs.select, 2, -1)
        self.assertRaises(ValueError, g.vs.select, (2, -1))
        self.assertRaises(ValueError, g.vs.__getitem__, (0, 1000000))
 
    def testGraphMethodProxying(self):
        g = Graph.Barabasi(100)
        vs = g.vs(1,3,5,7,9)
        self.assertEquals(vs.degree(), g.degree(vs))
        self.assertEquals(g.degree(vs), g.degree(vs.indices))
        for v, d in zip(vs, vs.degree()):
            self.assertEquals(v.degree(), d)


def suite():
    vertex_suite = unittest.makeSuite(VertexTests)
    vs_suite = unittest.makeSuite(VertexSeqTests)
    return unittest.TestSuite([vertex_suite, vs_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

