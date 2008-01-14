import unittest
from igraph import *

class MaxFlowTests(unittest.TestCase):
    def setUp(self):
        self.g=Graph(4, [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)])
        self.capacities = [4, 2, 10, 2, 2]
        for idx in range(self.g.ecount()):
            self.g.es[idx]["capacity"]=self.capacities[idx]

    def testCapacities(self):
        self.failUnless(self.capacities == \
            self.g.es.get_attribute_values("capacity"))

    def testMaxFlowValue(self):
        self.failUnless(self.g.maxflow_value(0, 3) == 2)
        self.failUnless(self.g.maxflow_value(0, 3, self.capacities) == 4)
        self.failUnless(self.g.maxflow_value(0, 3, "capacity") == 4)
        self.assertRaises(KeyError, self.g.maxflow_value, 0, 3, "unknown")

    def testMinCutValue(self):
        self.failUnless(self.g.mincut_value(0, 3) == 2)
        self.failUnless(self.g.mincut_value(0, 3, self.capacities) == 4)
        self.failUnless(self.g.mincut_value(0, 3, "capacity") == 4)
        self.assertRaises(KeyError, self.g.mincut_value, 0, 3, "unknown")
        self.failUnless(self.g.mincut_value() == 2)
        self.failUnless(self.g.mincut_value(source=0) == 2)
        self.failUnless(self.g.mincut_value(target=2) == 2)

def suite():
    flow_suite = unittest.makeSuite(MaxFlowTests)
    return unittest.TestSuite([flow_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

