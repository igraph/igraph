# vim:set ts=4 sw=4 sts=4 et:
import unittest
from igraph import *

class SpectralTests(unittest.TestCase):
    def testLaplacian(self):
        g=Graph.Full(3)
        g.es["weight"] = [1, 2, 3]
        self.failUnless(g.laplacian() == [[ 2, -1, -1],\
                                          [-1,  2, -1],\
                                          [-1, -1,  2]])
        self.failUnless(g.laplacian(normalized=True) == [[ 1.0, -0.5, -0.5],\
                                                        [-0.5,  1.0, -0.5],\
                                                        [-0.5, -0.5,  1.0]])
        self.failUnless(g.laplacian("weight") == [[ 3, -1, -2],\
                                                  [-1,  4, -3],\
                                                  [-2, -3,  5]])

        mx0 = [[1., -1/(12**0.5), -2/(15**0.5)],
              [-1/(12**0.5), 1., -3/(20**0.5)],
              [-2/(15**0.5), -3/(20**0.5), 1.]]
        mx1 = g.laplacian("weight", True)
        for i in xrange(3):
            for j in xrange(3):
                self.failUnless(abs(mx0[i][j]-mx1[i][j]) < 0.001)

        g=Graph.Tree(5, 2)
        g.add_vertices(1)
        self.failUnless(g.laplacian() == [[ 2, -1, -1,  0,  0, 0],\
                                          [-1,  3,  0, -1, -1, 0],\
                                          [-1,  0,  1,  0,  0, 0],\
                                          [ 0, -1,  0,  1,  0, 0],\
                                          [ 0, -1,  0,  0,  1, 0],\
                                          [ 0,  0,  0,  0,  0, 0]])
        
def suite():
    spectral_suite = unittest.makeSuite(SpectralTests)
    return unittest.TestSuite([spectral_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

