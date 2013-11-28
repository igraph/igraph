import random
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
    def tearDown(self):
        set_random_number_generator(random)

    def testSetRandomNumberGenerator(self):
        set_random_number_generator(FakeRNG)
        graph = Graph.GRG(10, 0.2)
        self.assertEqual(graph.vs["x"], [0.1] * 10)
        self.assertEqual(graph.vs["y"], [0.1] * 10)

        self.assertRaises(AttributeError, set_random_number_generator,
                InvalidRNG)

    def testSeeding(self):
        state = random.getstate()
        g1 = Graph.Erdos_Renyi(n=1000, m=5000)
        random.setstate(state)
        g2 = Graph.Erdos_Renyi(n=1000, m=5000)
        self.assertTrue(g1.get_edgelist() == g2.get_edgelist())


def suite():
    random_suite = unittest.makeSuite(RandomNumberGeneratorTests)
    return unittest.TestSuite([random_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

