import unittest
from igraph import *

class DirectedUndirectedTests(unittest.TestCase):
    def testToUndirected(self):
        graph = Graph([(0,1), (0,2), (1,0)], directed=True)

        graph2 = graph.copy()
        graph2.to_undirected(mode=False)
        self.assertTrue(graph2.vcount() == graph.vcount())
        self.assertTrue(graph2.is_directed() == False)
        self.assertTrue(sorted(graph2.get_edgelist()) == [(0,1), (0,1), (0,2)])

        graph2 = graph.copy()
        graph2.to_undirected()
        self.assertTrue(graph2.vcount() == graph.vcount())
        self.assertTrue(graph2.is_directed() == False)
        self.assertTrue(sorted(graph2.get_edgelist()) == [(0,1), (0,2)])

        graph2 = graph.copy()
        graph2.es["weight"] = [1,2,3]
        graph2.to_undirected(mode="collapse", combine_edges="sum")
        self.assertTrue(graph2.vcount() == graph.vcount())
        self.assertTrue(graph2.is_directed() == False)
        self.assertTrue(sorted(graph2.get_edgelist()) == [(0,1), (0,2)])
        self.assertTrue(graph2.es["weight"] == [4,2])

        graph = Graph([(0,1),(1,0),(0,1),(1,0),(2,1),(1,2)], directed=True)
        graph2 = graph.copy()
        graph2.es["weight"] = [1,2,3,4,5,6]
        graph2.to_undirected(mode="mutual", combine_edges="sum")
        self.assertTrue(graph2.vcount() == graph.vcount())
        self.assertTrue(graph2.is_directed() == False)
        self.assertTrue(sorted(graph2.get_edgelist()) == [(0,1), (0,1), (1,2)])
        self.assertTrue(graph2.es["weight"] == [7,3,11] or graph2.es["weight"] == [3,7,11])

    def testToDirected(self):
        graph = Graph([(0,1), (0,2), (2,3), (2,4)], directed=False)
        graph.to_directed()
        self.assertTrue(graph.is_directed())
        self.assertTrue(graph.vcount() == 5)
        self.assertTrue(sorted(graph.get_edgelist()) == \
                [(0,1), (0,2), (1,0), (2,0), (2,3), (2,4), (3,2), (4,2)]
        )


class GraphRepresentationTests(unittest.TestCase):
    def testGetAdjacency(self):
        # Undirected case
        g = Graph.Tree(6, 3)
        g.es["weight"] = range(5)
        self.assertTrue(g.get_adjacency() == Matrix([
            [0, 1, 1, 1, 0, 0],
            [1, 0, 0, 0, 1, 1],
            [1, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0]
        ]))
        self.assertTrue(g.get_adjacency(attribute="weight") == Matrix([
            [0, 0, 1, 2, 0, 0],
            [0, 0, 0, 0, 3, 4],
            [1, 0, 0, 0, 0, 0],
            [2, 0, 0, 0, 0, 0],
            [0, 3, 0, 0, 0, 0],
            [0, 4, 0, 0, 0, 0]
        ]))
        self.assertTrue(g.get_adjacency(eids=True) == Matrix([
            [0, 1, 2, 3, 0, 0],
            [1, 0, 0, 0, 4, 5],
            [2, 0, 0, 0, 0, 0],
            [3, 0, 0, 0, 0, 0],
            [0, 4, 0, 0, 0, 0],
            [0, 5, 0, 0, 0, 0]
        ])-1)

        # Directed case
        g = Graph.Tree(6, 3, "tree_out")
        g.add_edges([(0,1), (1,0)])
        self.assertTrue(g.get_adjacency() == Matrix([
            [0, 2, 1, 1, 0, 0],
            [1, 0, 0, 0, 1, 1],
            [0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0]
        ]))


def suite():
    direction_suite = unittest.makeSuite(DirectedUndirectedTests)
    representation_suite = unittest.makeSuite(GraphRepresentationTests)
    return unittest.TestSuite([direction_suite,
        representation_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

