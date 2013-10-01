import unittest

from igraph import *

def powerset(iterable):
    items_powers = [(item, 1 << i) for i, item in enumerate(iterable)]
    for i in range(1 << len(items_powers)):
        for item, power in items_powers:
            if i & power:
                yield item


leda_graph = Graph([
    (0,8),(0,12),(0,14),(1,9),(1,10),(1,13),
    (2,8),(2,9),(3,10),(3,11),(3,13),(4,9),(4,14),
    (5,14),(6,9),(6,14),(7,8),(7,12),(7,14)])
leda_graph.vs["type"] = [0]*8+[1]*7

class MatchingTests(unittest.TestCase):
    def setUp(self):
        self.matching = Matching(leda_graph,
                [12, 10, 8, 13, -1, 14, 9, -1, 2, 6, 1, -1, 0, 3, 5],
                "type")

    def testIsMaximal(self):
        self.assertTrue(self.matching.is_maximal())
        self.matching.matching[0] = -1
        self.matching.matching[12] = -1
        self.assertFalse(self.matching.is_maximal())

    def testMatchingRetrieval(self):
        m = [12, 10, 8, 13, -1, 14, 9, -1, 2, 6, 1, -1, 0, 3, 5]
        self.assertEqual(self.matching.matching, m)
        for i, mate in enumerate(m):
            if mate == -1:
                self.assertFalse(self.matching.is_matched(i))
                self.assertEqual(self.matching.match_of(i), None)
            else:
                self.assertTrue(self.matching.is_matched(i))
                self.assertEqual(self.matching.match_of(i), mate)
                self.assertEqual(self.matching.match_of(
                    leda_graph.vs[i]).index, leda_graph.vs[mate].index)


class MaximumBipartiteMatchingTests(unittest.TestCase):
    def testBipartiteMatchingSimple(self):
        # Specifying the "type" attribute explicitly
        matching = leda_graph.maximum_bipartite_matching("type")
        self.assertEqual(len(matching), 6)
        self.assertTrue(matching.is_maximal())

        # Using the default attribute
        matching = leda_graph.maximum_bipartite_matching()
        self.assertEqual(len(matching), 6)
        self.assertTrue(matching.is_maximal())

    def testBipartiteMatchingErrors(self):
        # Type vector too short
        g = Graph([(0, 1), (1, 2), (2, 3)])
        self.assertRaises(InternalError, g.maximum_bipartite_matching,
                types=[0,1,0])

        # Graph not bipartite
        self.assertRaises(ValueError, g.maximum_bipartite_matching,
                types=[0,1,1,1])


def suite():
    matching_suite = unittest.makeSuite(MatchingTests)
    bipartite_unweighted_suite = unittest.makeSuite(MaximumBipartiteMatchingTests)
    return unittest.TestSuite([matching_suite, bipartite_unweighted_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

