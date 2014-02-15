# vim:ts=4 sw=4 sts=4:
import unittest
from igraph import *

class AttributeTests(unittest.TestCase):
    def testGraphAttributes(self):
        g = Graph.Full(5)
        g["date"] = "2005-12-17"
        self.assertTrue(g["date"] == "2005-12-17")
        del g["date"]
        self.assertRaises(KeyError, g.__getitem__, "date")

    def testVertexAttributes(self):
        g = Graph.Full(5)
        g.vs[0]["name"] = "first"
        self.assertTrue(g.vs[0]["name"] == "first")
        del g.vs["name"]
        self.assertRaises(KeyError, g.vs.__getitem__, "name")

        g.vs[0]["name"] = "second"
        g.vs[0]["date"] = "2007-12-17"
        ans = g.vs[0].attribute_names()
        ans.sort()
        self.assertTrue(ans == ["date", "name"])
        attrs = g.vs[0].attributes()
        self.assertTrue(attrs == {"name": "second", "date": "2007-12-17"})

    def testEdgeAttributes(self):
        g = Graph.Full(5)
        g.es[0]["name"] = "first"
        self.assertTrue(g.es[0]["name"] == "first")
        del g.es["name"]
        self.assertRaises(KeyError, g.es.__getitem__, "name")

        g.es[0]["name"] = "second"
        g.es[0]["date"] = "2007-12-17"
        ans = g.es[0].attribute_names()
        ans.sort()
        self.assertTrue(ans == ["date", "name"])
        attrs = g.es[0].attributes()
        self.assertTrue(attrs == {"name": "second", "date": "2007-12-17"})

    def testMassVertexAttributeAssignment(self):
        g = Graph.Full(5)
        g.vs.set_attribute_values("name", range(5))
        self.assertTrue(g.vs.get_attribute_values("name") == range(5))
        g.vs["name"] = range(5,10)
        self.assertTrue(g.vs["name"] == range(5,10))
        g.vs["name2"] = (1,2,3,4,6) 
        self.assertTrue(g.vs["name2"] == [1,2,3,4,6])
        g.vs.set_attribute_values("name", [2])
        self.assertTrue(g.vs["name"] == [2]*5)

    def testMassEdgeAttributeAssignment(self):
        g = Graph.Full(5)
        g.es.set_attribute_values("name", range(10))
        self.assertTrue(g.es.get_attribute_values("name") == range(10))
        g.es["name"] = range(10,20)
        self.assertTrue(g.es["name"] == range(10,20))
        g.es["name2"] = (1,2,3,4,6,1,2,3,4,6) 
        self.assertTrue(g.es["name2"] == [1,2,3,4,6,1,2,3,4,6])
        g.es.set_attribute_values("name", [2])
        self.assertTrue(g.es["name"] == [2]*10)

    def testVertexNameIndexing(self):
        g = Graph.Famous("bull")
        g.vs["name"] = ["foo", "bar", "baz", "fred", "thud"]
        self.assertTrue(g.degree("bar") == 3)
        self.assertTrue(g.degree(["bar", "fred", 0]) == [3, 1, 2])
        g.vs[2]["name"] = "quack"
        self.assertRaises(ValueError, g.degree, "baz")
        self.assertTrue(g.degree("quack") == 3)
        self.assertTrue(g.degree(u"quack") == 3)
        self.assertTrue(g.degree([u"bar", u"thud", 0]) == [3, 1, 2])
        del g.vs["name"]
        self.assertRaises(ValueError, g.degree, [u"bar", u"thud", 0])

    def testInvalidAttributeNames(self):
        g = Graph.Famous("bull")
        for attr_name in [None, 2.654, unittest, str]:
            self.assertRaises(TypeError, g.vs.__setitem__, attr_name, "foo")
            self.assertRaises(TypeError, g.vs.__getitem__, attr_name, "foo")
            self.assertRaises(TypeError, g.vs[0].__setitem__, attr_name, "foo")
            self.assertRaises(TypeError, g.vs[0].__getitem__, attr_name, "foo")
            self.assertRaises(TypeError, g.es.__setitem__, attr_name, "foo")
            self.assertRaises(TypeError, g.es.__getitem__, attr_name, "foo")
            self.assertRaises(TypeError, g.es[0].__setitem__, attr_name, "foo")
            self.assertRaises(TypeError, g.es[0].__getitem__, attr_name, "foo")

class AttributeCombinationTests(unittest.TestCase):
    def setUp(self):
        el = [(0,1), (1,0), (1,2), (2,3), (2,3), (2,3), (3,3)] 
        self.g = Graph(el)
        self.g.es["weight"] = [1, 2, 3, 4, 5, 6, 7]
        self.g.es["weight2"] = [1, 2, 3, 4, 5, 6, 7]

    def testCombinationMax(self):
        g = self.g
        g.simplify(combine_edges="max")
        self.assertTrue(g.es["weight"] == [2, 3, 6])
        self.assertTrue(g.es["weight2"] == [2, 3, 6])

    def testCombinationMin(self):
        g = self.g
        g.simplify(combine_edges="min")
        self.assertTrue(g.es["weight"] == [1, 3, 4])
        self.assertTrue(g.es["weight2"] == [1, 3, 4])

    def testCombinationRandom(self):
        g = self.g
        g.simplify(combine_edges="random")
        del g.es["weight2"]
        for i in xrange(100):
            self.assertTrue(g.es[0]["weight"] in (1, 2))
            self.assertTrue(g.es[1]["weight"] == 3)
            self.assertTrue(g.es[2]["weight"] in (4, 5, 6))

    def testCombinationMean(self):
        g = self.g
        g.simplify(combine_edges="mean")
        self.assertTrue(g.es["weight"] == [1.5, 3, 5])
        self.assertTrue(g.es["weight2"] == [1.5, 3, 5])

    def testCombinationMedian(self):
        g = self.g
        g.es["weight2"] = [1, 0, 2, 4, 8, 6, 7]
        g.simplify(combine_edges="median")
        self.assertTrue(g.es["weight"] == [1.5, 3, 5])
        self.assertTrue(g.es["weight2"] == [0.5, 2, 6])

    def testCombinationSum(self):
        g = self.g
        g.simplify(combine_edges="sum")
        self.assertTrue(g.es["weight"] == [3, 3, 15])
        self.assertTrue(g.es["weight2"] == [3, 3, 15])

    def testCombinationProd(self):
        g = self.g
        g.simplify(combine_edges="prod")
        self.assertTrue(g.es["weight"] == [2, 3, 120])
        self.assertTrue(g.es["weight2"] == [2, 3, 120])

    def testCombinationMedian(self):
        g = self.g
        g.es["weight2"] = [1, 0, 2, 4, 8, 6, 7]
        g.simplify(combine_edges="median")
        self.assertTrue(g.es["weight"] == [1.5, 3, 5])
        self.assertTrue(g.es["weight2"] == [0.5, 2, 6])

    def testCombinationFirst(self):
        g = self.g
        g.es["weight2"] = [1, 0, 2, 6, 8, 4, 7]
        g.simplify(combine_edges="first")
        self.assertTrue(g.es["weight"] == [1, 3, 4])
        self.assertTrue(g.es["weight2"] == [1, 2, 6])

    def testCombinationLast(self):
        g = self.g
        g.es["weight2"] = [1, 0, 2, 6, 8, 4, 7]
        g.simplify(combine_edges="last")
        self.assertTrue(g.es["weight"] == [2, 3, 6])
        self.assertTrue(g.es["weight2"] == [0, 2, 4])

    def testCombinationConcat(self):
        g = self.g
        g.es["name"] = list("ABCDEFG")
        g.simplify(combine_edges=dict(name="concat"))
        self.assertFalse("weight" in g.edge_attributes())
        self.assertFalse("weight2" in g.edge_attributes())
        self.assertTrue(g.es["name"] == ["AB", "C", "DEF"])

    def testCombinationMaxMinIgnore(self):
        g = self.g
        g.es["name"] = list("ABCDEFG")
        g.simplify(combine_edges={"weight": "min", "weight2": "max", "name": "ignore"})
        self.assertTrue(g.es["weight"] == [1, 3, 4])
        self.assertTrue(g.es["weight2"] == [2, 3, 6])
        self.assertFalse("name" in g.edge_attributes())

    def testCombinationIgnoreAsNone(self):
        g = self.g
        g.es["name"] = list("ABCDEFG")
        g.simplify(combine_edges={"weight": "min", "name": None})
        self.assertTrue(g.es["weight"] == [1, 3, 4])
        self.assertFalse("weight2" in g.edge_attributes())
        self.assertFalse("name" in g.edge_attributes())

    def testCombinationFunction(self):
        g = self.g

        def join_dash(l):
            return "-".join(l)

        g.es["name"] = list("ABCDEFG")
        g.simplify(combine_edges={"weight": max, "name": join_dash})
        self.assertTrue(g.es["weight"] == [2, 3, 6])
        self.assertFalse("weight2" in g.edge_attributes())
        self.assertTrue(g.es["name"] == ["A-B", "C", "D-E-F"])

    def testCombinationDefault(self):
        g = self.g
        g.simplify(combine_edges={None: "max"})
        self.assertTrue(g.es["weight"] == [2, 3, 6])
        self.assertTrue(g.es["weight2"] == [2, 3, 6])

    def testCombinationDefaultOverride(self):
        g = self.g
        g.simplify(combine_edges={None: "max", "weight": "sum"})
        self.assertTrue(g.es["weight"] == [3, 3, 15])
        self.assertTrue(g.es["weight2"] == [2, 3, 6])

    def testCombinationTypeMismatch(self):
        g = self.g
        g.es["weight"] = list("ABCDEFG")
        self.assertRaises(TypeError, g.simplify, combine_edges={"weight": "mean"})

    def testCombinationNonexistentAttribute(self):
        g = self.g
        g.simplify(combine_edges={"nonexistent": max})
        self.assertTrue(g.edge_attributes() == [])

    def testCombinationNone(self):
        g = self.g
        g.simplify()
        self.assertTrue(sorted(g.edge_attributes()) == [])


class UnicodeAttributeTests(unittest.TestCase):
    def testUnicodeAttributeNameCombination(self):
        g = Graph.Erdos_Renyi(n=9, m=20)
        g.es[u"test"] = 1
        g.contract_vertices([0,0,0,1,1,1,2,2,2])


def suite():
    attribute_suite = unittest.makeSuite(AttributeTests)
    attribute_combination_suite = unittest.makeSuite(AttributeCombinationTests)
    unicode_attributes_suite = unittest.makeSuite(UnicodeAttributeTests)
    return unittest.TestSuite([attribute_suite, attribute_combination_suite,
        unicode_attributes_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()
