import unittest
from igraph import Graph, Layout

class LayoutTests(unittest.TestCase):
    def testConstructor(self):
        layout = Layout([(0,0,1), (0,1,0), (1,0,0)])
        self.assertEqual(layout.dim, 3)
        layout = Layout([(0,0,1), (0,1,0), (1,0,0)], 3)
        self.assertEqual(layout.dim, 3)
        self.assertRaises(ValueError, Layout, [(0,1), (1,0)], 3)
        self.assertRaises(ValueError, Layout, [])

    def testIndexing(self):
        layout = Layout([(0,0,1), (0,1,0), (1,0,0), (2,1,3)])
        self.assertEqual(len(layout), 4)
        self.assertEqual(layout[1], [0, 1, 0])
        self.assertEqual(layout[3], [2, 1, 3])
        
        row = layout[2]
        row[2] = 1
        self.assertEqual(layout[2], [1, 0, 0])
        layout[2] = row
        self.assertEqual(layout[2], [1, 0, 1])

        del layout[1]
        self.assertEqual(len(layout), 3)

    def testScaling(self):
        layout = Layout([(0,0,1), (0,1,0), (1,0,0), (2,1,3)])
        layout.scale(1.5)
        self.assertEqual(layout.coords, [[0., 0., 1.5], \
                                         [0., 1.5, 0.], \
                                         [1.5, 0., 0.], \
                                         [3., 1.5, 4.5]])

        layout = Layout([(0,0,1), (0,1,0), (1,0,0), (2,1,3)])
        layout.scale(1, 1, 3)
        self.assertEqual(layout.coords, [[0, 0, 3], \
                                         [0, 1, 0], \
                                         [1, 0, 0], \
                                         [2, 1, 9]])


        layout = Layout([(0,0,1), (0,1,0), (1,0,0), (2,1,3)])
        layout.scale((2, 2, 1))
        self.assertEqual(layout.coords, [[0, 0, 1], \
                                         [0, 2, 0], \
                                         [2, 0, 0], \
                                         [4, 2, 3]])

        self.assertRaises(ValueError, layout.scale, 2, 3)


    def testTranslation(self):
        layout = Layout([(0,0,1), (0,1,0), (1,0,0), (2,1,3)])
        layout2 = layout.copy()
        
        layout.translate(1,3,2)
        self.assertEqual(layout.coords, [[1, 3, 3], \
                                         [1, 4, 2], \
                                         [2, 3, 2], \
                                         [3, 4, 5]])
        
        layout.translate((-1,-3,-2))
        self.assertEqual(layout.coords, layout2.coords)

        self.assertRaises(ValueError, layout.translate, v=[3])
        
    def testCentroid(self):
        layout = Layout([(0,0,1), (0,1,0), (1,0,0), (2,1,3)])
        centroid = layout.centroid()
        self.assertEqual(len(centroid), 3)
        self.assertAlmostEqual(centroid[0], 0.75)
        self.assertAlmostEqual(centroid[1], 0.5)
        self.assertAlmostEqual(centroid[2], 1.)


    def testBoundingBox(self):
        layout = Layout([(0,0,1), (0,1,0), (1,0,0), (2,1,3)])
        self.assertEqual(layout.bounding_box(), (0,0,0,2,1,3))
        self.assertEqual(layout.bounding_box(1), (-1,-1,-1,3,2,4))


    def testCenter(self):
        layout = Layout([(-2,0), (-2,-2), (0,-2), (0,0)])
        layout.center()
        self.assertEqual(layout.coords, [[-1,1], [-1,-1], [1,-1], [1,1]])
        layout.center(5,5)
        self.assertEqual(layout.coords, [[4,6], [4,4], [6,4], [6,6]])
        self.assertRaises(ValueError, layout.center, 3)
        self.assertRaises(TypeError, layout.center, p=6)


    def testToPolar(self):
        import math
        layout = Layout([(0, 0), (-1, 1), (0, 1), (1, 1)])
        layout.to_radial(min_angle = 180, max_angle = 0, max_radius = 2)
        exp = [[0., 0.], [-2., 0.], [0., 2.], [2, 0.]]
        for idx in xrange(4):
            self.assertAlmostEqual(layout.coords[idx][0], exp[idx][0], places=3)
            self.assertAlmostEqual(layout.coords[idx][1], exp[idx][1], places=3)

    def testTransform(self):
        def tr(coord, dx, dy): return coord[0]+dx, coord[1]+dy
        layout = Layout([(1, 2), (3, 4)])
        layout.transform(tr, 2, -1)
        self.assertEqual(layout.coords, [[3, 1], [5, 3]])


class LayoutAlgorithmTests(unittest.TestCase):
    def testFruchtermanReingold(self):
        g = Graph.Barabasi(100)
        lo = g.layout("fr")
        self.failUnless(isinstance(lo, Layout))
        lo = g.layout("fr", miny=range(100))
        self.failUnless(isinstance(lo, Layout))
        lo = g.layout("fr", miny=range(100), maxy=range(100))
        self.failUnless(isinstance(lo, Layout))

    def testReingoldTilford(self):
        g = Graph.Barabasi(100)
        lo = g.layout("rt")
        ys = [coord[1] for coord in lo]
        root = ys.index(0.0)
        self.assertEqual(ys, g.shortest_paths(root)[0])
        g = Graph.Barabasi(100) + Graph.Barabasi(50)
        lo = g.layout("rt", root=[0, 100])
        self.assertEqual(lo[100][1]-lo[0][1], 0)
        lo = g.layout("rt", root=[0, 100], rootlevel = [2, 10])
        self.assertEqual(lo[100][1]-lo[0][1], 8)


def suite():
    layout_suite = unittest.makeSuite(LayoutTests)
    layout_algorithm_suite = unittest.makeSuite(LayoutAlgorithmTests)
    return unittest.TestSuite([layout_suite, layout_algorithm_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

