import unittest
from igraph import *

class IteratorTests(unittest.TestCase):
    def testBFS(self):
        g=Graph.Tree(10, 2)
        vs=[v.index for v in g.bfsiter(0)]
        self.assertEqual(vs, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        vs=[(v.index,dist,parent) for v,dist,parent in g.bfsiter(0, advanced=True)]
        vs=[(v,d,p.index) for v,d,p in vs if p != None]
        self.assertEqual(vs, [(1,1,0), (2,1,0), (3,2,1), (4,2,1), \
          (5,2,2), (6,2,2), (7,3,3), (8,3,3), (9,3,4)])


def suite():
    iterator_suite = unittest.makeSuite(IteratorTests)
    return unittest.TestSuite([iterator_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

