import unittest

from itertools import izip
from igraph import *

class ColorTests(unittest.TestCase):
    def assertAlmostEqualMany(self, items1, items2, eps):
        for idx, (item1, item2) in enumerate(izip(items1, items2)):
            self.assertAlmostEqual(item1, item2, places=eps,
                    msg="mismatch at index %d, %r != %r with %d digits"
                    % (idx, items1, items2, eps))

    def setUp(self):
        columns = ["r", "g", "b", "h", "v", "l", "s_hsv", "s_hsl", "alpha"]
        # Examples taken from http://en.wikipedia.org/wiki/HSL_and_HSV
        values = [
                (1, 1, 1, 0, 1, 1, 0, 0, 1),
                (0.5, 0.5, 0.5, 0, 0.5, 0.5, 0, 0, 0.5),
                (0, 0, 0, 0, 0, 0, 0, 0, 1),
                (1, 0, 0, 0, 1, 0.5, 1, 1, 0.5),
                (0.75, 0.75, 0, 60, 0.75, 0.375, 1, 1, 0.25),
                (0, 0.5, 0, 120, 0.5, 0.25, 1, 1, 0.75),
                (0.5, 1, 1, 180, 1, 0.75, 0.5, 1, 1),
                (0.5, 0.5, 1, 240, 1, 0.75, 0.5, 1, 1),
                (0.75, 0.25, 0.75, 300, 0.75, 0.5, 0.666666667, 0.5, 0.25),
                (0.211, 0.149, 0.597, 248.3, 0.597, 0.373, 0.750, 0.601, 1),
                (0.495, 0.493, 0.721, 240.5, 0.721, 0.607, 0.316, 0.290, 0.75),
        ]
        self.data = [dict(zip(columns, value)) for value in values]
        for row in self.data:
            row["h"] /= 360.

    def _testGeneric(self, method, args1, args2=("r", "g", "b")):
        if len(args1) == len(args2)+1:
            args2 += ("alpha", )
        for data in self.data:
            vals1 = [data.get(arg, 0.0) for arg in args1]
            vals2 = [data.get(arg, 0.0) for arg in args2]
            self.assertAlmostEqualMany(method(*vals1), vals2, 2)

    def testHSVtoRGB(self):
        self._testGeneric(hsv_to_rgb, "h s_hsv v".split())

    def testHSVAtoRGBA(self):
        self._testGeneric(hsva_to_rgba, "h s_hsv v alpha".split())

    def testHSLtoRGB(self):
        self._testGeneric(hsl_to_rgb, "h s_hsl l".split())

    def testHSLAtoRGBA(self):
        self._testGeneric(hsla_to_rgba, "h s_hsl l alpha".split())

    def testRGBtoHSL(self):
        self._testGeneric(rgb_to_hsl, "r g b".split(), "h s_hsl l".split())

    def testRGBAtoHSLA(self):
        self._testGeneric(rgba_to_hsla, "r g b alpha".split(), "h s_hsl l alpha".split())

    def testRGBtoHSV(self):
        self._testGeneric(rgb_to_hsv, "r g b".split(), "h s_hsv v".split())

    def testRGBAtoHSVA(self):
        self._testGeneric(rgba_to_hsva, "r g b alpha".split(), "h s_hsv v alpha".split())


class PaletteTests(unittest.TestCase):
    def testGradientPalette(self):
        gp = GradientPalette("red", "blue", 3)
        self.assertTrue(gp.get(0) == (1., 0., 0., 1.))
        self.assertTrue(gp.get(1) == (0.5, 0., 0.5, 1.))
        self.assertTrue(gp.get(2) == (0., 0., 1., 1.))

    def testAdvancedGradientPalette(self):
        agp = AdvancedGradientPalette(["red", "black", "blue"], n=9)
        self.assertTrue(agp.get(0) == (1., 0., 0., 1.))
        self.assertTrue(agp.get(2) == (0.5, 0., 0., 1.))
        self.assertTrue(agp.get(4) == (0., 0., 0., 1.))
        self.assertTrue(agp.get(5) == (0., 0., 0.25, 1.))
        self.assertTrue(agp.get(8) == (0., 0., 1., 1.))

        agp = AdvancedGradientPalette(["red", "black", "blue"], [0, 8, 2], 9)
        self.assertTrue(agp.get(0) == (1., 0., 0., 1.))
        self.assertTrue(agp.get(1) == (0.5, 0., 0.5, 1.))
        self.assertTrue(agp.get(5) == (0., 0., 0.5, 1.))


def suite():
    color_suite = unittest.makeSuite(ColorTests)
    palette_suite = unittest.makeSuite(PaletteTests)
    return unittest.TestSuite([color_suite, palette_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

