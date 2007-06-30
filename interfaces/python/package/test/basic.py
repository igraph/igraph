import unittest
from igraph import *

class BasicTests(unittest.TestCase):
    def testGraphCreation(self):
        self.failUnless(isinstance(Graph(), Graph))
        g=Graph(3, [(0,1), (1,2), (2,0)])
        self.failUnless(g.vcount() == 3 and g.ecount() == 3 and g.is_directed() == False)
        g=Graph(2, [(0,1), (1,2), (2,3)], True)
        self.failUnless(g.vcount() == 4 and g.ecount() == 3 and g.is_directed() == True)
        g=Graph([(0,1), (1,2)])
        self.failUnless(g.vcount() == 3 and g.ecount() == 2 and g.is_directed() == False)

    def testPickling(self):
        import pickle
        g=Graph([(0,1), (1,2)])
        g["data"]="abcdef"
        g.vs["data"]=[3,4,5]
        g.es["data"]=["A", "B"]
        pickled=pickle.dumps(g)
        self.failUnless(isinstance(pickled, str))
        g2=pickle.loads(pickled)
        self.failUnless(g["data"] == g2["data"])
        self.failUnless(g.vs["data"] == g2.vs["data"])
        self.failUnless(g.es["data"] == g2.es["data"])
        self.failUnless(g.vcount() == g2.vcount())
        self.failUnless(g.ecount() == g2.ecount())
        self.failUnless(g.is_directed() == g2.is_directed())


def suite():
    basic_suite = unittest.makeSuite(BasicTests)
    return unittest.TestSuite([basic_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

