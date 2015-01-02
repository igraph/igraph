import random
import unittest
from igraph import Graph, InternalError


class RandomWalkTests(unittest.TestCase):
    def validate_walk(self, g, walk, start, length, mode="out"):
        prev = None
        for vertex in walk:
            if prev is not None:
                self.assertTrue(vertex in g.neighbors(prev, mode=mode))
            else:
                self.assertEquals(start, vertex)
            prev = vertex

    def testRandomWalkUndirected(self):
        g = Graph.GRG(100, 0.2)
        for i in xrange(100):
            start = random.randint(0, g.vcount()-1)
            length = random.randint(0, 10)
            walk = g.random_walk(start, length)
            self.validate_walk(g, walk, start, length)

    def testRandomWalkDirectedOut(self):
        g = Graph.Tree(121, 3, mode="out")
        mode = "out"
        for i in xrange(100):
            start = 0
            length = random.randint(0, 4)
            walk = g.random_walk(start, length, mode)
            self.validate_walk(g, walk, start, length, mode)

    def testRandomWalkDirectedIn(self):
        g = Graph.Tree(121, 3, mode="out")
        mode = "in"
        for i in xrange(100):
            start = random.randint(40, g.vcount()-1)
            length = random.randint(0, 4)
            walk = g.random_walk(start, length, mode)
            self.validate_walk(g, walk, start, length, mode)

    def testRandomWalkDirectedAll(self):
        g = Graph.Tree(121, 3, mode="out")
        mode = "all"
        for i in xrange(100):
            start = random.randint(0, g.vcount()-1)
            length = random.randint(0, 10)
            walk = g.random_walk(start, length, mode)
            self.validate_walk(g, walk, start, length, mode)

    def testRandomWalkStuck(self):
        g = Graph.Ring(10, circular=False, directed=True)
        walk = g.random_walk(5, 20)
        self.assertEquals([5, 6, 7, 8, 9], walk)
        self.assertRaises(InternalError, g.random_walk, 5, 20, stuck="error")


def suite():
    random_walk_suite = unittest.makeSuite(RandomWalkTests)
    return unittest.TestSuite([random_walk_suite])


def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())


if __name__ == "__main__":
    test()
