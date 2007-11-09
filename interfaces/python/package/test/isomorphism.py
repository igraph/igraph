import unittest
from igraph import *

class IsomorphismTests(unittest.TestCase):
    def testIsomorphic(self):
        g1 = Graph(8, [(0, 4), (0, 5), (0, 6), \
                       (1, 4), (1, 5), (1, 7), \
                       (2, 4), (2, 6), (2, 7), \
                       (3, 5), (3, 6), (3, 7)])
        g2 = Graph(8, [(0, 1), (0, 3), (0, 4), \
                       (2, 3), (2, 1), (2, 6), \
                       (5, 1), (5, 4), (5, 6), \
                       (7, 3), (7, 6), (7, 4)])
        self.failUnless(g1.isomorphic(g2))
        self.failUnless(g2.isomorphic(g1, return_mapping_21=True) == (True, None, [0, 2, 5, 7, 1, 3, 4, 6]))

def suite():
    isomorphism_suite = unittest.makeSuite(IsomorphismTests)
    return unittest.TestSuite([isomorphism_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

