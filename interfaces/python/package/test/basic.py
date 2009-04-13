import unittest
from igraph import *

class BasicTests(unittest.TestCase):
    def testGraphCreation(self):
        self.failUnless(isinstance(Graph(), Graph))
        g=Graph(3, [(0,1), (1,2), (2,0)])
        self.failUnless(g.vcount() == 3 and g.ecount() == 3 and g.is_directed() == False and g.is_simple() == True)
        g=Graph(2, [(0,1), (1,2), (2,3)], True)
        self.failUnless(g.vcount() == 4 and g.ecount() == 3 and g.is_directed() == True and g.is_simple())
        g=Graph([(0,1), (1,2), (2,1)])
        self.failUnless(g.vcount() == 3 and g.ecount() == 3 and g.is_directed() == False and g.is_simple() == False)
        g=Graph([(0,1), (0,0), (1,2)])
        self.failUnless(g.vcount() == 3 and g.ecount() == 3 and g.is_directed() == False and g.is_simple() == False)


    def testAdjacency(self):
        g=Graph(4, [(0,1), (1,2), (2,0), (2,3)], directed=True)
        self.failUnless(g.neighbors(2) == [0, 1, 3])
        self.failUnless(g.predecessors(2) == [1])
        self.failUnless(g.successors(2) == [0, 3])
        self.failUnless(g.get_adjlist() == [[1], [2], [0,3], []])
        self.failUnless(g.get_adjlist(IN) == [[2], [0], [1], [2]])
        self.failUnless(g.get_adjlist(ALL) == [[1,2], [0,2], [0,1,3], [2]])
        
    def testEdgeAdjacency(self):
        g=Graph(4, [(0,1), (1,2), (2,0), (2,3)], directed=True)
        self.failUnless(g.adjacent(2) == [2, 3])
        self.failUnless(g.adjacent(2, IN) == [1])
        self.failUnless(g.adjacent(2, ALL) == [2, 3, 1])
        self.failUnless(g.get_adjedgelist() == [[0], [1], [2,3], []])
        self.failUnless(g.get_adjedgelist(IN) == [[2], [0], [1], [3]])
        self.failUnless(g.get_adjedgelist(ALL) == [[0,2], [1,0], [2,3,1], [3]])


    def testMultiplesLoops(self):
        g=Graph.Tree(7, 2)
        g.add_vertices(1)
        g.add_edges([(0,1), (7,7), (6,6), (6,6), (6,6)])
        
        # is_loop
        self.failUnless(g.is_loop() == [False, False, False, False, \
            False, False, False, True, True, True, True])
        self.failUnless(g.is_loop(g.ecount()-2))
        self.failUnless(g.is_loop(range(6,8)) == [False, True])

        # is_multiple
        self.failUnless(g.is_multiple() == [False, False, False, False, \
            False, False, True, False, False, True, True])

        # count_multiple
        self.failUnless(g.count_multiple() == [2, 1, 1, 1, 1, 1, 2, 1, 3, 3, 3])
        self.failUnless(g.count_multiple(g.ecount()-1) == 3)
        self.failUnless(g.count_multiple(range(2,5)) == [1, 1, 1])

        # check if a mutual directed edge pair is reported as multiple
        g=Graph(2, [(0,1), (1,0)], directed=True)
        self.failUnless(g.is_multiple() == [False, False])


    def testPickling(self):
        import pickle
        g=Graph([(0,1), (1,2)])
        g["data"]="abcdef"
        g.vs["data"]=[3,4,5]
        g.es["data"]=["A", "B"]
        pickled=pickle.dumps(g)
        self.failUnless(isinstance(pickled, str))
        g2=pickle.loads(pickled)
        self.failUnless(g["data"] == g2["data"])
        self.failUnless(g.vs["data"] == g2.vs["data"])
        self.failUnless(g.es["data"] == g2.es["data"])
        self.failUnless(g.vcount() == g2.vcount())
        self.failUnless(g.ecount() == g2.ecount())
        self.failUnless(g.is_directed() == g2.is_directed())


class DatatypeTests(unittest.TestCase):
    def testMatrix(self):
        m = Matrix([[1,2,3], [4,5], [6,7,8]])
        self.failUnless(m.shape == (3, 3))

        # Reading data
        self.failUnless(m.data == [[1,2,3], [4,5,0], [6,7,8]])
        self.failUnless(m[1,1] == 5)
        self.failUnless(m[0] == [1,2,3])
        self.failUnless(m[0,:] == [1,2,3])
        self.failUnless(m[:,0] == [1,4,6])
        self.failUnless(m[2,0:2] == [6,7])
        self.failUnless(m[:,:].data == [[1,2,3], [4,5,0], [6,7,8]])
        self.failUnless(m[:,1:3].data == [[2,3], [5,0], [7,8]])

        # Writing data
        m[1,1] = 10
        self.failUnless(m[1,1] == 10)
        m[1] = (6,5,4)
        self.failUnless(m[1] == [6,5,4])
        m[1:3] = [[4,5,6], (7,8,9)]
        self.failUnless(m[1:3].data == [[4,5,6], [7,8,9]])

        # Minimums and maximums
        self.failUnless(m.min() == 1)
        self.failUnless(m.max() == 9)
        self.failUnless(m.min(0) == [1,2,3])
        self.failUnless(m.max(0) == [7,8,9])
        self.failUnless(m.min(1) == [1,4,7])
        self.failUnless(m.max(1) == [3,6,9])

        # Special constructors
        m=Matrix.Fill(2, (3,3))
        self.failUnless(m.min() == 2 and m.max() == 2 and m.shape == (3,3))
        m=Matrix.Zero(5, 4)
        self.failUnless(m.min() == 0 and m.max() == 0 and m.shape == (5,4))
        m=Matrix.Identity(3)
        self.failUnless(m.data == [[1,0,0], [0,1,0], [0,0,1]])
        m=Matrix.Identity(3, 2)
        self.failUnless(m.data == [[1,0], [0,1], [0,0]])


class GraphDictListTests(unittest.TestCase):
    def setUp(self):
        self.vertices = [ \
                {"name": "Alice", "age": 48, "gender": "F"},
                {"name": "Bob",   "age": 33, "gender": "M"},
                {"name": "Cecil", "age": 45, "gender": "F"},
                {"name": "David", "age": 34, "gender": "M"}
        ]
        self.edges = [ \
                {"source": "Alice", "target": "Bob",   "friendship": 4, "advice": 4},
                {"source": "Cecil", "target": "Bob",   "friendship": 5, "advice": 5},
                {"source": "Cecil", "target": "Alice", "friendship": 5, "advice": 5},
                {"source": "David", "target": "Alice", "friendship": 2, "advice": 4},
                {"source": "David", "target": "Bob",   "friendship": 1, "advice": 2}
        ]

    def testGraphFromDictList(self):
        g = Graph.DictList(self.vertices, self.edges)
        self.checkIfOK(g, "name")
        g = Graph.DictList(self.vertices, self.edges, iterative=True)
        self.checkIfOK(g, "name")

    def testGraphFromDictIterator(self):
        g = Graph.DictList(iter(self.vertices), iter(self.edges))
        self.checkIfOK(g, "name")
        g = Graph.DictList(iter(self.vertices), iter(self.edges), iterative=True)
        self.checkIfOK(g, "name")

    def testGraphFromDictIteratorNoVertices(self):
        g = Graph.DictList(None, iter(self.edges))
        self.checkIfOK(g, "name", check_vertex_attrs=False)
        g = Graph.DictList(None, iter(self.edges), iterative=True)
        self.checkIfOK(g, "name", check_vertex_attrs=False)

    def testGraphFromDictListExtraVertexName(self):
        del self.vertices[2:]      # No data for "Cecil" and "David"
        g = Graph.DictList(self.vertices, self.edges)
        self.failUnless(g.vcount() == 4 and g.ecount() == 5 and g.is_directed() == False)
        self.failUnless(g.vs["name"] == ["Alice", "Bob", "Cecil", "David"])
        self.failUnless(g.vs["age"] == [48, 33, None, None])
        self.failUnless(g.vs["gender"] == ["F", "M", None, None])
        self.failUnless(g.es["friendship"] == [4, 5, 5, 2, 1])
        self.failUnless(g.es["advice"] == [4, 5, 5, 4, 2])
        self.failUnless(g.get_edgelist() == [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3)])

    def testGraphFromDictListAlternativeName(self):
        for vdata in self.vertices:
            vdata["name_alternative"] = vdata["name"]
            del vdata["name"]
        g = Graph.DictList(self.vertices, self.edges, vertex_name_attr="name_alternative")
        self.checkIfOK(g, "name_alternative")
        g = Graph.DictList(self.vertices, self.edges, vertex_name_attr="name_alternative", \
                iterative=True)
        self.checkIfOK(g, "name_alternative")

    def checkIfOK(self, g, name_attr, check_vertex_attrs=True):
        self.failUnless(g.vcount() == 4 and g.ecount() == 5 and g.is_directed() == False)
        self.failUnless(g.get_edgelist() == [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3)])
        self.failUnless(g.vs[name_attr] == ["Alice", "Bob", "Cecil", "David"])
        if check_vertex_attrs:
            self.failUnless(g.vs["age"] == [48, 33, 45, 34])
            self.failUnless(g.vs["gender"] == ["F", "M", "F", "M"])
        self.failUnless(g.es["friendship"] == [4, 5, 5, 2, 1])
        self.failUnless(g.es["advice"] == [4, 5, 5, 4, 2])


def suite():
    basic_suite = unittest.makeSuite(BasicTests)
    datatype_suite = unittest.makeSuite(DatatypeTests)
    graph_dict_list_suite = unittest.makeSuite(GraphDictListTests)
    return unittest.TestSuite([basic_suite, datatype_suite, graph_dict_list_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

