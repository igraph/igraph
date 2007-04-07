import unittest
from igraph import *

class DecompositionTests(unittest.TestCase):
    def testKCores(self):
        g = Graph(11, [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3),
                       (2,4), (2,5), (3,6), (3,7), (1,7), (7,8),
                       (1,9), (1,10), (9,10)])
        self.failUnless(g.coreness() == [3,3,3,3,1,1,1,2,1,2,2])
        self.failUnless(g.shell_index() == g.coreness())

        l=g.k_core(3).get_edgelist()
        l.sort()
        self.failUnless(l == [(0,1), (0,2), (0,3), (1,2), (1,3), (2,3)])

class CommunityTests(unittest.TestCase):
    def testClauset(self):
        g = Graph.Full(5) + Graph.Full(5)
        g.add_edges([(0, 5)])
        cl, q = g.community_clauset(return_q = True)
        self.failUnless(cl == [0,0,0,0,0,1,1,1,1,1])
        self.assertAlmostEqual(q, 0.4523)
        self.failUnless(max(g.community_clauset(3)) == 2)
        
def suite():
    decomposition_suite = unittest.makeSuite(DecompositionTests)
    return unittest.TestSuite([decomposition_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

