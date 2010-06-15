import unittest
from igraph import *


class FakeRNG(object):
    @staticmethod
    def random():
        return 0.1

    @staticmethod
    def randint(a, b):
        return a

    @staticmethod
    def gauss(mu, sigma):
        return 0.3

class InvalidRNG(object):
    pass


class RandomNumberGeneratorTests(unittest.TestCase):
    def testSetRandomNumberGenerator(self):
        set_random_number_generator(FakeRNG)
        graph, xs, ys = Graph.GRG(10, 0.2, return_coordinates=True)
        self.failUnless(xs == ys and xs == [0.1] * 10)

        self.assertRaises(AttributeError, set_random_number_generator,
                InvalidRNG)

def suite():
    random_suite = unittest.makeSuite(RandomNumberGeneratorTests)
    return unittest.TestSuite([random_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

