import unittest
from igraph import *

class BasicTests(unittest.TestCase):
    def testGraphCreation(self):
        g = Graph()
        self.failUnless(isinstance(g, Graph))
        self.failUnless(g.vcount() == 0 and g.ecount() == 0 and g.is_directed() == False)

        g=Graph(3, [(0,1), (1,2), (2,0)])
        self.failUnless(g.vcount() == 3 and g.ecount() == 3 and g.is_directed() == False and g.is_simple() == True)
        g=Graph(2, [(0,1), (1,2), (2,3)], True)
        self.failUnless(g.vcount() == 4 and g.ecount() == 3 and g.is_directed() == True and g.is_simple())
        g=Graph([(0,1), (1,2), (2,1)])
        self.failUnless(g.vcount() == 3 and g.ecount() == 3 and g.is_directed() == False and g.is_simple() == False)
        g=Graph([(0,1), (0,0), (1,2)])
        self.failUnless(g.vcount() == 3 and g.ecount() == 3 and g.is_directed() == False and g.is_simple() == False)

    def testAddVertex(self):
        g = Graph()
        g.add_vertex()
        self.failUnless(g.vcount() == 1 and g.ecount() == 0)
        self.failIf("name" in g.vertex_attributes())
        g.add_vertex("foo")
        self.failUnless(g.vcount() == 2 and g.ecount() == 0)
        self.failUnless("name" in g.vertex_attributes())
        self.assertEquals(g.vs["name"], [None, "foo"])
        g.add_vertex(3)
        self.failUnless(g.vcount() == 3 and g.ecount() == 0)
        self.failUnless("name" in g.vertex_attributes())
        self.assertEquals(g.vs["name"], [None, "foo", 3])
        g.add_vertex(name="bar")
        self.failUnless(g.vcount() == 4 and g.ecount() == 0)
        self.failUnless("name" in g.vertex_attributes())
        self.assertEquals(g.vs["name"], [None, "foo", 3, "bar"])
        g.add_vertex(name="frob", spam="cheese", ham=42)
        self.failUnless(g.vcount() == 5 and g.ecount() == 0)
        self.assertEquals(sorted(g.vertex_attributes()), ["ham", "name", "spam"])
        self.assertEquals(g.vs["spam"], [None]*4 + ["cheese"])
        self.assertEquals(g.vs["ham"], [None]*4 + [42])

    def testAddVertices(self):
        g = Graph()
        g.add_vertices(2)
        self.failUnless(g.vcount() == 2 and g.ecount() == 0)
        g.add_vertices("spam")
        self.failUnless(g.vcount() == 3 and g.ecount() == 0)
        self.assertEquals(g.vs[2]["name"], "spam")
        g.add_vertices(["bacon", "eggs"])
        self.failUnless(g.vcount() == 5 and g.ecount() == 0)
        self.assertEquals(g.vs[2:]["name"], ["spam", "bacon", "eggs"])

    def testDeleteVertices(self):
        g = Graph([(0,1), (1,2), (2,3), (0,2), (3,4), (4,5)])
        self.assertEquals(6, g.vcount())
        self.assertEquals(6, g.ecount())

        # Delete a single vertex
        g.delete_vertices(4)
        self.assertEquals(5, g.vcount())
        self.assertEquals(4, g.ecount())

        # Delete multiple vertices
        g.delete_vertices([1,3])
        self.assertEquals(3, g.vcount())
        self.assertEquals(1, g.ecount())

        # Delete a vertex sequence
        g.delete_vertices(g.vs[:2])
        self.assertEquals(1, g.vcount())
        self.assertEquals(0, g.ecount())

        # Delete a single vertex object
        g.vs[0].delete()
        self.assertEquals(0, g.vcount())
        self.assertEquals(0, g.ecount())

        # Delete vertices by name
        g = Graph.Full(4)
        g.vs["name"] = ["spam", "bacon", "eggs", "ham"]
        self.assertEquals(4, g.vcount())
        g.delete_vertices("spam")
        self.assertEquals(3, g.vcount())
        g.delete_vertices(["bacon", "ham"])
        self.assertEquals(1, g.vcount())

        # Deleting a nonexistent vertex
        self.assertRaises(ValueError, g.delete_vertices, "no-such-vertex")
        self.assertRaises(InternalError, g.delete_vertices, 2)

    def testAddEdges(self):
        g = Graph()
        g.add_vertices(["spam", "bacon", "eggs", "ham"])

        g.add_edges([(0, 1)])
        self.assertEquals(g.vcount(), 4)
        self.assertEquals(g.get_edgelist(), [(0, 1)])

        g.add_edges([(1, 2), (2, 3), (1, 3)])
        self.assertEquals(g.vcount(), 4)
        self.assertEquals(g.get_edgelist(), [(0, 1), (1, 2), (2, 3), (1, 3)])

        g.add_edges([("spam", "eggs"), ("spam", "ham")])
        self.assertEquals(g.vcount(), 4)
        self.assertEquals(g.get_edgelist(), [(0, 1), (1, 2), (2, 3), (1, 3), (0, 2), (0, 3)])

    def testDeleteEdges(self):
        g = Graph.Famous("petersen")
        g.vs["name"] = list("ABCDEFGHIJ")
        el = g.get_edgelist()

        self.assertEquals(15, g.ecount())

        # Deleting single edge
        g.delete_edges(14)
        el[14:] = []
        self.assertEquals(14, g.ecount())
        self.assertEquals(el, g.get_edgelist())

        # Deleting multiple edges
        g.delete_edges([2,5,7])
        el[7:8] = []; el[5:6] = []; el[2:3] = []
        self.assertEquals(11, g.ecount())
        self.assertEquals(el, g.get_edgelist())

        # Deleting edge object
        g.es[6].delete()
        el[6:7] = []
        self.assertEquals(10, g.ecount())
        self.assertEquals(el, g.get_edgelist())

        # Deleting edge sequence object
        g.es[1:4].delete()
        el[1:4] = []
        self.assertEquals(7, g.ecount())
        self.assertEquals(el, g.get_edgelist())

        # Deleting edges by IDs
        g.delete_edges([(2,7), (5,8)])
        el[4:5] = []; el[1:2] = []
        self.assertEquals(5, g.ecount())
        self.assertEquals(el, g.get_edgelist())

        # Deleting edges by names
        g.delete_edges([("D", "I"), ("G", "I")])
        el[3:4] = []; el[1:2] = []
        self.assertEquals(3, g.ecount())
        self.assertEquals(el, g.get_edgelist())

        # Deleting nonexistent edges
        self.assertRaises(ValueError, g.delete_edges, [(0,2)])
        self.assertRaises(ValueError, g.delete_edges, [("A", "C")])
        self.assertRaises(ValueError, g.delete_edges, [(0,15)])

    def testGraphGetEids(self):
        g = Graph.Famous("petersen")
        eids = g.get_eids(pairs=[(0,1), (0,5), (1, 6), (4, 9), (8, 6)])
        self.failUnless(eids == [0, 2, 4, 9, 12])
        eids = g.get_eids(path=[0, 1, 2, 3, 4])
        self.failUnless(eids == [0, 3, 5, 7])
        eids = g.get_eids(pairs=[(7,9), (9,6)], path = [7, 9, 6])
        self.failUnless(eids == [14, 13, 14, 13])
        self.assertRaises(InternalError, g.get_eids, pairs=[(0,1), (0,2)])

    def testAdjacency(self):
        g=Graph(4, [(0,1), (1,2), (2,0), (2,3)], directed=True)
        self.failUnless(g.neighbors(2) == [0, 1, 3])
        self.failUnless(g.predecessors(2) == [1])
        self.failUnless(g.successors(2) == [0, 3])
        self.failUnless(g.get_adjlist() == [[1], [2], [0,3], []])
        self.failUnless(g.get_adjlist(IN) == [[2], [0], [1], [2]])
        self.failUnless(g.get_adjlist(ALL) == [[1,2], [0,2], [0,1,3], [2]])
        
    def testEdgeIncidency(self):
        g=Graph(4, [(0,1), (1,2), (2,0), (2,3)], directed=True)
        self.failUnless(g.incident(2) == [2, 3])
        self.failUnless(g.incident(2, IN) == [1])
        self.failUnless(g.incident(2, ALL) == [2, 3, 1])
        self.failUnless(g.get_inclist() == [[0], [1], [2,3], []])
        self.failUnless(g.get_inclist(IN) == [[2], [0], [1], [3]])
        self.failUnless(g.get_inclist(ALL) == [[0,2], [1,0], [2,3,1], [3]])


    def testMultiplesLoops(self):
        g=Graph.Tree(7, 2)

        # has_multiple
        self.failIf(g.has_multiple())

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

        # has_multiple
        self.failUnless(g.has_multiple())

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
        g.custom_data = None
        pickled=pickle.dumps(g)

        g2=pickle.loads(pickled)
        self.failUnless(g["data"] == g2["data"])
        self.failUnless(g.vs["data"] == g2.vs["data"])
        self.failUnless(g.es["data"] == g2.es["data"])
        self.failUnless(g.vcount() == g2.vcount())
        self.failUnless(g.ecount() == g2.ecount())
        self.failUnless(g.is_directed() == g2.is_directed())
        self.failUnless(g2.custom_data == g.custom_data)


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

        # Conversion to string
        m=Matrix.Identity(3)
        self.failUnless(str(m) == "[[1, 0, 0]\n [0, 1, 0]\n [0, 0, 1]]")
        self.failUnless(repr(m) == "Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])")


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


class GraphTupleListTests(unittest.TestCase):
    def setUp(self):
        self.edges = [ \
                ("Alice", "Bob", 4, 4),
                ("Cecil", "Bob", 5, 5),
                ("Cecil", "Alice", 5, 5),
                ("David", "Alice", 2, 4),
                ("David", "Bob", 1, 2)
        ]

    def testGraphFromTupleList(self):
        g = Graph.TupleList(self.edges)
        self.checkIfOK(g, "name", ())

    def testGraphFromTupleListWithEdgeAttributes(self):
        g = Graph.TupleList(self.edges, edge_attrs=("friendship", "advice"))
        self.checkIfOK(g, "name", ("friendship", "advice"))
        g = Graph.TupleList(self.edges, edge_attrs=("friendship", ))
        self.checkIfOK(g, "name", ("friendship", ))
        g = Graph.TupleList(self.edges, edge_attrs="friendship")
        self.checkIfOK(g, "name", ("friendship", ))

    def testGraphFromTupleListWithDifferentNameAttribute(self):
        g = Graph.TupleList(self.edges, vertex_name_attr="spam")
        self.checkIfOK(g, "spam", ())

    def testGraphFromTupleListWithWeights(self):
        g = Graph.TupleList(self.edges, weights=True)
        self.checkIfOK(g, "name", ("weight", ))
        g = Graph.TupleList(self.edges, weights="friendship")
        self.checkIfOK(g, "name", ("friendship", ))
        g = Graph.TupleList(self.edges, weights=False)
        self.checkIfOK(g, "name", ())
        self.assertRaises(ValueError, Graph.TupleList,
                [self.edges], weights=True, edge_attrs="friendship")

    def testNoneForMissingAttributes(self):
        g = Graph.TupleList(self.edges, edge_attrs=("friendship", "advice", "spam"))
        self.checkIfOK(g, "name", ("friendship", "advice", "spam"))

    def checkIfOK(self, g, name_attr, edge_attrs):
        self.failUnless(g.vcount() == 4 and g.ecount() == 5 and g.is_directed() == False)
        self.failUnless(g.get_edgelist() == [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3)])
        self.failUnless(g.attributes() == [])
        self.failUnless(g.vertex_attributes() == [name_attr])
        self.failUnless(g.vs[name_attr] == ["Alice", "Bob", "Cecil", "David"])
        if edge_attrs:
            self.failUnless(sorted(g.edge_attributes()) == sorted(edge_attrs))
            self.failUnless(g.es[edge_attrs[0]] == [4, 5, 5, 2, 1])
            if len(edge_attrs) > 1:
                self.failUnless(g.es[edge_attrs[1]] == [4, 5, 5, 4, 2])
            if len(edge_attrs) > 2:
                self.failUnless(g.es[edge_attrs[2]] == [None] * 5)
        else:
            self.failUnless(g.edge_attributes() == [])


class DegreeSequenceTests(unittest.TestCase):
    def testIsDegreeSequence(self):
        self.assertTrue(is_degree_sequence([]))
        self.assertTrue(is_degree_sequence([], []))
        self.assertTrue(is_degree_sequence([0]))
        self.assertTrue(is_degree_sequence([0], [0]))
        self.assertFalse(is_degree_sequence([1]))
        self.assertTrue(is_degree_sequence([1], [1]))
        self.assertTrue(is_degree_sequence([2]))
        self.assertFalse(is_degree_sequence([2, 1, 1, 1]))
        self.assertTrue(is_degree_sequence([2, 1, 1, 1], [1, 1, 1, 2]))
        self.assertFalse(is_degree_sequence([2, 1, -2]))
        self.assertFalse(is_degree_sequence([2, 1, 1, 1], [1, 1, 1, -2]))
        self.assertTrue(is_degree_sequence([3, 3, 3, 3, 3, 3, 3, 3, 3, 3]))
        self.assertTrue(is_degree_sequence([3, 3, 3, 3, 3, 3, 3, 3, 3, 3], None))
        self.assertFalse(is_degree_sequence([3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]))
        self.assertTrue(is_degree_sequence([3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
            [4, 3, 2, 3, 4, 4, 2, 2, 4, 2]))

    def testIsGraphicalSequence(self):
        self.assertTrue(is_graphical_degree_sequence([]))
        self.assertTrue(is_graphical_degree_sequence([], []))
        self.assertTrue(is_graphical_degree_sequence([0]))
        self.assertTrue(is_graphical_degree_sequence([0], [0]))
        self.assertFalse(is_graphical_degree_sequence([1]))
        self.assertTrue(is_graphical_degree_sequence([1], [1]))
        self.assertFalse(is_graphical_degree_sequence([2]))
        self.assertFalse(is_graphical_degree_sequence([2, 1, 1, 1]))
        self.assertTrue(is_graphical_degree_sequence([2, 1, 1, 1], [1, 1, 1, 2]))
        self.assertFalse(is_graphical_degree_sequence([2, 1, -2]))
        self.assertFalse(is_graphical_degree_sequence([2, 1, 1, 1], [1, 1, 1, -2]))
        self.assertTrue(is_graphical_degree_sequence([3, 3, 3, 3, 3, 3, 3, 3, 3, 3]))
        self.assertTrue(is_graphical_degree_sequence([3, 3, 3, 3, 3, 3, 3, 3, 3, 3], None))
        self.assertFalse(is_graphical_degree_sequence([3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]))
        self.assertTrue(is_graphical_degree_sequence([3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
            [4, 3, 2, 3, 4, 4, 2, 2, 4, 2]))
        self.assertTrue(is_graphical_degree_sequence([3, 3, 3, 3, 4]))


def suite():
    basic_suite = unittest.makeSuite(BasicTests)
    datatype_suite = unittest.makeSuite(DatatypeTests)
    graph_dict_list_suite = unittest.makeSuite(GraphDictListTests)
    graph_tuple_list_suite = unittest.makeSuite(GraphTupleListTests)
    degree_sequence_suite = unittest.makeSuite(DegreeSequenceTests)
    return unittest.TestSuite([basic_suite, datatype_suite, graph_dict_list_suite,
        graph_tuple_list_suite, degree_sequence_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

