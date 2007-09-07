import unittest
from igraph import *

class GeneratorTests(unittest.TestCase):
    def testFamous(self):
        g=Graph.Famous("tutte")
        self.failUnless(g.vcount() == 46 and g.ecount() == 69)
        self.assertRaises(InternalError, Graph.Famous, "unknown")

    def testFull(self):
        g=Graph.Full(20, directed=True)
        el=g.get_edgelist()
        el.sort()
        self.failUnless(g.get_edgelist() == [(x, y) for x in range(20) for y in range(20) if x!=y])

    def testLCF(self):
        g1=Graph.LCF(12, (5, -5), 6)
        g2=Graph.Famous("Franklin")
        self.failUnless(g1.isomorphic(g2))
        self.assertRaises(InternalError, Graph.LCF, 12, (5, -5), -3)

    def testKautz(self):
        g=Graph.Kautz(2, 2)
        deg_in=g.degree(type=IN)
        deg_out=g.degree(type=OUT)
        # This is not a proper test, but should spot most errors
        self.failUnless(g.is_directed() and deg_in==[2]*12 and deg_out==[2]*12)

    def testDeBruijn(self):
        g=Graph.De_Bruijn(2, 3)
        deg_in=g.degree(type=IN, loops=True)
        deg_out=g.degree(type=OUT, loops=True)
        # This is not a proper test, but should spot most errors
        self.failUnless(g.is_directed() and deg_in==[2]*8 and deg_out==[2]*8)

        
def suite():
    generator_suite = unittest.makeSuite(GeneratorTests)
    return unittest.TestSuite([generator_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

