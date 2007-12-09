import unittest
from igraph import *

class ColorTests(unittest.TestCase):
    def testGradientPalette(self):
        gp = GradientPalette("red", "blue", 3)
        self.failUnless(gp.get(0) == (1., 0., 0.))
        self.failUnless(gp.get(1) == (0.5, 0., 0.5))
        self.failUnless(gp.get(2) == (0., 0., 1.))

    def testAdvancedGradientPalette(self):
        agp = AdvancedGradientPalette(["red", "black", "blue"], n=9)
        self.failUnless(agp.get(0) == (1., 0., 0.))
        self.failUnless(agp.get(2) == (0.5, 0., 0.))
        self.failUnless(agp.get(4) == (0., 0., 0.))
        self.failUnless(agp.get(5) == (0., 0., 0.25))
        self.failUnless(agp.get(8) == (0., 0., 1.))

        agp = AdvancedGradientPalette(["red", "black", "blue"], [0, 8, 2], 9)
        self.failUnless(agp.get(0) == (1., 0., 0.))
        self.failUnless(agp.get(1) == (0.5, 0., 0.5))
        self.failUnless(agp.get(5) == (0., 0., 0.5))


def suite():
    color_suite = unittest.makeSuite(ColorTests)
    return unittest.TestSuite([color_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

