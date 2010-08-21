import unittest
from igraph import *

class SpectralTests(unittest.TestCase):
    def assertAlmostEqualsMatrix(self, mat1, mat2, eps = 1e-7):
        self.failUnless(all(
            abs(obs-exp) < eps
            for obs, exp in zip(sum(mat1, []), sum(mat2, []))
        ))

    def testLaplacian(self):
        g=Graph.Full(3)
        self.failUnless(g.laplacian() == [[ 2, -1, -1],\
                                          [-1,  2, -1],\
                                          [-1, -1,  2]])
        self.assertAlmostEqualsMatrix(g.laplacian(normalized=True),
                [[1, -0.5, -0.5], [-0.5, 1, -0.5], [-0.5, -0.5, 1]])

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

