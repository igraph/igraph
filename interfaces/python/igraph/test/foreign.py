import unittest
from igraph import *
import tempfile
import os

class ForeignTests(unittest.TestCase):
    def testDIMACS(self):
        tmpf, tmpfname = tempfile.mkstemp()
        os.close(tmpf)
        tmpf = open(tmpfname, "w")
        print >>tmpf, """c
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
"""
        tmpf.close()
        g, src, dst, cap=Graph.Read_DIMACS(tmpfname, False)
        self.failUnless(isinstance(g, Graph))
        self.failUnless(g.vcount() == 4 and g.ecount() == 5)
        self.failUnless(src == 0 and dst == 3)
        self.failUnless(cap == [4,2,2,3,5])

        g.write_dimacs(tmpfname, src, dst, cap)
        os.unlink(tmpfname)


    def testAdjacency(self):
        tmpf, tmpfname = tempfile.mkstemp()
        os.close(tmpf)
        tmpf = open(tmpfname, "w")
        print >>tmpf, """# Test comment line
0 1 1 0 0 0
1 0 1 0 0 0
1 1 0 0 0 0
0 0 0 0 2 2
0 0 0 2 0 2
0 0 0 2 2 0
"""
        tmpf.close()
        g = Graph.Read_Adjacency(tmpfname)
        self.failUnless(isinstance(g, Graph))
        self.failUnless(g.vcount() == 6 and g.ecount() == 18 and
            g.is_directed() and "weight" not in g.edge_attributes())
        g = Graph.Read_Adjacency(tmpfname, attribute="weight")
        self.failUnless(isinstance(g, Graph))
        self.failUnless(g.vcount() == 6 and g.ecount() == 12 and
            g.is_directed() and g.es["weight"] == [1,1,1,1,1,1,2,2,2,2,2,2])

        g.write_adjacency(tmpfname)
        os.unlink(tmpfname)


def suite():
    foreign_suite = unittest.makeSuite(ForeignTests)
    return unittest.TestSuite([foreign_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

