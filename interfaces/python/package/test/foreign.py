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
p min 4 5
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
	

def suite():
    foreign_suite = unittest.makeSuite(ForeignTests)
    return unittest.TestSuite([foreign_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

