# vim:ts=4 sw=4 sts=4:
import unittest
from igraph import *

class AttributeTests(unittest.TestCase):
    def testGraphAttributes(self):
        g = Graph.Full(5)
        g["date"] = "2005-12-17"
        self.failUnless(g["date"] == "2005-12-17")
        del g["date"]
        self.assertRaises(KeyError, g.__getitem__, "date")

    def testVertexAttributes(self):
        g = Graph.Full(5)
        g.vs[0]["name"] = "first"
        self.failUnless(g.vs[0]["name"] == "first")
        del g.vs["name"]
        self.assertRaises(KeyError, g.vs.__getitem__, "name")

        g.vs[0]["name"] = "second"
        g.vs[0]["date"] = "2007-12-17"
        ans = g.vs[0].attribute_names()
        ans.sort()
        self.failUnless(ans == ["date", "name"])
        attrs = g.vs[0].attributes()
        self.failUnless(attrs == {"name": "second", "date": "2007-12-17"})

    def testEdgeAttributes(self):
        g = Graph.Full(5)
        g.es[0]["name"] = "first"
        self.failUnless(g.es[0]["name"] == "first")
        del g.es["name"]
        self.assertRaises(KeyError, g.es.__getitem__, "name")

        g.es[0]["name"] = "second"
        g.es[0]["date"] = "2007-12-17"
        ans = g.es[0].attribute_names()
        ans.sort()
        self.failUnless(ans == ["date", "name"])
        attrs = g.es[0].attributes()
        self.failUnless(attrs == {"name": "second", "date": "2007-12-17"})

    def testMassVertexAttributeAssignment(self):
        g = Graph.Full(5)
        g.vs.set_attribute_values("name", range(5))
        self.failUnless(g.vs.get_attribute_values("name") == range(5))
        g.vs["name"] = range(5,10)
        self.failUnless(g.vs["name"] == range(5,10))
        g.vs["name2"] = (1,2,3,4,6) 
        self.failUnless(g.vs["name2"] == [1,2,3,4,6])
        g.vs.set_attribute_values("name", [2])
        self.failUnless(g.vs["name"] == [2]*5)

    def testMassEdgeAttributeAssignment(self):
        g = Graph.Full(5)
        g.es.set_attribute_values("name", range(10))
        self.failUnless(g.es.get_attribute_values("name") == range(10))
        g.es["name"] = range(10,20)
        self.failUnless(g.es["name"] == range(10,20))
        g.es["name2"] = (1,2,3,4,6,1,2,3,4,6) 
        self.failUnless(g.es["name2"] == [1,2,3,4,6,1,2,3,4,6])
        g.es.set_attribute_values("name", [2])
        self.failUnless(g.es["name"] == [2]*10)

    def testVertexNameIndexing(self):
        g = Graph.Famous("bull")
        g.vs["name"] = ["foo", "bar", "baz", "fred", "thud"]
        self.failUnless(g.degree("bar") == 3)
        self.failUnless(g.degree(["bar", "fred", 0]) == [3, 1, 2])
        g.vs[2]["name"] = "quack"
        self.assertRaises(ValueError, g.degree, "baz")
        self.failUnless(g.degree("quack") == 3)
        self.failUnless(g.degree(u"quack") == 3)
        self.failUnless(g.degree([u"bar", u"thud", 0]) == [3, 1, 2])
        del g.vs["name"]
        self.assertRaises(ValueError, g.degree, [u"bar", u"thud", 0])

    
def suite():
    attribute_suite = unittest.makeSuite(AttributeTests)
    return unittest.TestSuite([attribute_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()
