from __future__ import with_statement

import unittest
from igraph import *
from igraph.test.utils import temporary_file


class ForeignTests(unittest.TestCase):
    def testDIMACS(self):
        with temporary_file(u"""\
        c
        c        This is a simple example file to demonstrate the
        c     DIMACS input file format for minimum-cost flow problems.
        c
        c problem line :
        p max 4 5
        c
        c node descriptor lines :
        n 1 s
        n 4 t
        c
        c arc descriptor lines :
        a 1 2 4
        a 1 3 2
        a 2 3 2
        a 2 4 3
        a 3 4 5
        """) as tmpfname:
            graph = Graph.Read_DIMACS(tmpfname, False)
            self.assertTrue(isinstance(graph, Graph))
            self.assertTrue(graph.vcount() == 4 and graph.ecount() == 5)
            self.assertTrue(graph["source"] == 0 and graph["target"] == 3)
            self.assertTrue(graph.es["capacity"] == [4,2,2,3,5])
            graph.write_dimacs(tmpfname)


    def testDL(self):
        with temporary_file(u"""\
        dl n=5
        format = fullmatrix
        labels embedded
        data:
        larry david lin pat russ
        Larry 0 1 1 1 0
        david 1 0 0 0 1
        Lin 1 0 0 1 0
        Pat 1 0 1 0 1
        russ 0 1 0 1 0
        """) as tmpfname:
            g = Graph.Read_DL(tmpfname)
            self.assertTrue(isinstance(g, Graph))
            self.assertTrue(g.vcount() == 5 and g.ecount() == 12)
            self.assertTrue(g.is_directed())
            self.assertTrue(sorted(g.get_edgelist()) == \
                    [(0,1),(0,2),(0,3),(1,0),(1,4),(2,0),(2,3),(3,0),\
                     (3,2),(3,4),(4,1),(4,3)])

        with temporary_file(u"""\
        dl n=5
        format = fullmatrix
        labels:
        barry,david
        lin,pat
        russ
        data:
        0 1 1 1 0
        1 0 0 0 1
        1 0 0 1 0
        1 0 1 0 1
        0 1 0 1 0
        """) as tmpfname:
            g = Graph.Read_DL(tmpfname)
            self.assertTrue(isinstance(g, Graph))
            self.assertTrue(g.vcount() == 5 and g.ecount() == 12)
            self.assertTrue(g.is_directed())
            self.assertTrue(sorted(g.get_edgelist()) == \
                    [(0,1),(0,2),(0,3),(1,0),(1,4),(2,0),(2,3),(3,0),\
                     (3,2),(3,4),(4,1),(4,3)])

        with temporary_file(u"""\
        DL n=5
        format = edgelist1
        labels:
        george, sally, jim, billy, jane
        labels embedded:
        data:
        george sally 2
        george jim 3
        sally jim 4
        billy george 5
        jane jim 6
        """) as tmpfname:
            g = Graph.Read_DL(tmpfname, False)
            self.assertTrue(isinstance(g, Graph))
            self.assertTrue(g.vcount() == 5 and g.ecount() == 5)
            self.assertTrue(not g.is_directed())
            self.assertTrue(sorted(g.get_edgelist()) == \
                    [(0,1),(0,2),(0,3),(1,2),(2,4)])

    def _testNCOLOrLGL(self, func, fname):
        g = func(fname, names=False, weights=False, \
                directed=False)
        self.assertTrue(isinstance(g, Graph))
        self.assertTrue(g.vcount() == 4 and g.ecount() == 5)
        self.assertTrue(not g.is_directed())
        self.assertTrue(sorted(g.get_edgelist()) == \
                [(0,1),(0,2),(1,1),(1,3),(2,3)])
        self.assertTrue("name" not in g.vertex_attributes() and \
                "weight" not in g.edge_attributes())

        g = func(fname, names=False, \
                directed=False)
        self.assertTrue("name" not in g.vertex_attributes() and \
                "weight" in g.edge_attributes())
        self.assertTrue(g.es["weight"] == [1, 2, 0, 3, 0])

        g = func(fname, directed=False)
        self.assertTrue("name" in g.vertex_attributes() and \
                "weight" in g.edge_attributes())
        self.assertTrue(g.vs["name"] == ["eggs", "spam", "ham", "bacon"])
        self.assertTrue(g.es["weight"] == [1, 2, 0, 3, 0])

    def testNCOL(self):
        with temporary_file(u"""\
        eggs spam 1
        ham eggs 2
        ham bacon
        bacon spam 3
        spam spam""") as tmpfname:
            self._testNCOLOrLGL(func=Graph.Read_Ncol, fname=tmpfname)

        with temporary_file(u"""\
        eggs spam
        ham eggs
        ham bacon
        bacon spam
        spam spam""") as tmpfname:
            g = Graph.Read_Ncol(tmpfname)
            self.assertTrue("name" in g.vertex_attributes() and \
                "weight" not in g.edge_attributes())

    def testLGL(self):
        with temporary_file(u"""\
        # eggs
        spam 1
        # ham
        eggs 2
        bacon
        # bacon
        spam 3
        # spam
        spam""") as tmpfname:
            self._testNCOLOrLGL(func=Graph.Read_Lgl, fname=tmpfname)

        with temporary_file(u"""\
        # eggs
        spam
        # ham
        eggs
        bacon
        # bacon
        spam
        # spam
        spam""") as tmpfname:
            g = Graph.Read_Lgl(tmpfname)
            self.assertTrue("name" in g.vertex_attributes() and \
                "weight" not in g.edge_attributes())


    def testAdjacency(self):
        with temporary_file(u"""\
        # Test comment line
        0 1 1 0 0 0
        1 0 1 0 0 0
        1 1 0 0 0 0
        0 0 0 0 2 2
        0 0 0 2 0 2
        0 0 0 2 2 0
        """) as tmpfname:
            g = Graph.Read_Adjacency(tmpfname)
            self.assertTrue(isinstance(g, Graph))
            self.assertTrue(g.vcount() == 6 and g.ecount() == 18 and
                g.is_directed() and "weight" not in g.edge_attributes())
            g = Graph.Read_Adjacency(tmpfname, attribute="weight")
            self.assertTrue(isinstance(g, Graph))
            self.assertTrue(g.vcount() == 6 and g.ecount() == 12 and
                g.is_directed() and g.es["weight"] == [1,1,1,1,1,1,2,2,2,2,2,2])

            g.write_adjacency(tmpfname)

    def testPickle(self):
        pickle = [128, 2, 99, 105, 103, 114, 97, 112, 104, 10, 71, 114, 97, 112,
                  104, 10, 113, 1, 40, 75, 3, 93, 113, 2, 75, 1, 75, 2, 134, 113, 3, 97,
                  137, 125, 125, 125, 116, 82, 113, 4, 125, 98, 46]
        if sys.version_info > (3, 0):
            pickle = bytes(pickle)
        else:
            pickle = "".join(map(chr, pickle))
        with temporary_file(pickle, "wb") as tmpfname:
            g = Graph.Read_Pickle(pickle)
            self.assertTrue(isinstance(g, Graph))
            self.assertTrue(g.vcount() == 3 and g.ecount() == 1 and
                not g.is_directed())
            g.write_pickle(tmpfname)


def suite():
    foreign_suite = unittest.makeSuite(ForeignTests)
    return unittest.TestSuite([foreign_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())

if __name__ == "__main__":
    test()

