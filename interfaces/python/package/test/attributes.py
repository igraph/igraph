import unittest
from igraph import *

class AttributeTests(unittest.TestCase):
    def testMassVertexAttributeAssignment(self):
	g = Graph.Full(5)
	g.vs.set_attribute_values("name", range(5))
	self.failUnless(g.vs.get_attribute_values("name") == range(5))
	self.assertRaises(ValueError, g.vs.set_attribute_values, "name", [2])

    def testMassEdgeAttributeAssignment(self):
	g = Graph.Full(5)
	g.es.set_attribute_values("name", range(10))
	self.failUnless(g.es.get_attribute_values("name") == range(10))
	self.assertRaises(ValueError, g.es.set_attribute_values, "name", [2])
	
def suite():
    attribute_suite = unittest.makeSuite(AttributeTests)
    return unittest.TestSuite([attribute_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()
