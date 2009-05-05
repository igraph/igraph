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

    def testEdgeConnectivity(self):
        self.failUnless(self.g.edge_connectivity(0, 3) == 2)
        self.failUnless(Graph.Barabasi(50).edge_connectivity() == 1)
        self.failUnless(self.g.adhesion() == 2)
        self.failUnless(Graph.Tree(10, 3).adhesion() == 1)
        self.failUnless(Graph.Tree(10, 3, TREE_OUT).adhesion() == 0)
        self.assertRaises(ValueError, self.g.edge_connectivity, 0)

    def testVertexConnectivity(self):
        self.failUnless(self.g.vertex_connectivity(0, 3) == 2)
        self.failUnless(Graph.Barabasi(50).vertex_connectivity() == 1)
        self.failUnless(self.g.cohesion() == 2)
        self.failUnless(Graph.Tree(10, 3).cohesion() == 1)
        self.failUnless(Graph.Tree(10, 3, TREE_OUT).cohesion() == 0)
        self.assertRaises(ValueError, self.g.vertex_connectivity, 0)
        self.assertRaises(InternalError, self.g.vertex_connectivity, 0, 1)
        self.failUnless(self.g.vertex_connectivity(0, 1, neighbors="infinity") > 1000)

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

    def testMinCut(self):
        mc = self.g.mincut()
        self.failUnless(isinstance(mc, Cut))
        self.failUnless(mc.value == 2)
        self.failUnless(set(mc.partition[0]).union(mc.partition[1]) == \
          set(range(self.g.vcount())))
        self.failUnless(isinstance(str(mc), str))
        self.failUnless(isinstance(repr(mc), str))
        self.failUnless(isinstance(mc.es, EdgeSeq))
        self.failUnless(len(mc.es) == 2)
        mc = self.g.mincut(self.capacities)
        self.failUnless(mc.value == 4)
        self.assertRaises(KeyError, self.g.mincut, "unknown")


def suite():
    flow_suite = unittest.makeSuite(MaxFlowTests)
    return unittest.TestSuite([flow_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

