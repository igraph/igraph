import unittest

from igraph import *

def powerset(iterable):
    items_powers = [(item, 1 << i) for i, item in enumerate(iterable)]
    for i in range(1 << len(items_powers)):
        for item, power in items_powers:
            if i & power:
                yield item

class IsSeparatorTests(unittest.TestCase):
    def testIsSeparator(self):
        g = Graph.Lattice([8, 4], circular=False)
        self.failUnless(g.is_separator([3, 11, 19, 27]))
        self.failIf(g.is_separator([10, 11, 18, 19]))
        self.failUnless(g.is_separator([29, 20, 11, 2]))
        self.failUnless(g.is_separator([16, 25, 17]))

        g = Graph.Lattice([8, 4], circular=True)
        self.failIf(g.is_separator([3, 11, 19, 27]))
        self.failIf(g.is_separator([29, 20, 11, 2]))

        self.assertRaises(InternalError, g.is_separator, range(32))

    def testIsMinimalSeparator(self):
        g = Graph.Lattice([8, 4], circular=False)
        self.failUnless(g.is_minimal_separator([3, 11, 19, 27]))
        self.failIf(g.is_minimal_separator([3, 11, 19, 27, 28]))
        self.failIf(g.is_minimal_separator([16, 25, 17]))
        self.failUnless(g.is_minimal_separator([16, 25]))

        self.assertRaises(InternalError, g.is_minimal_separator, range(32))

    def testAllMinimalSTSeparators(self):
        g = Graph.Famous("petersen")
        min_st_seps = set(tuple(x) for x in g.all_minimal_st_separators())
        for vs in powerset(range(g.vcount())):
            if vs in min_st_seps:
                self.failUnless(g.is_minimal_separator(vs))
            else:
                self.failIf(g.is_minimal_separator(vs))

    def testMinimumSizeSeparators(self):
        g = Graph.Famous("zachary")
        min_st_seps = set(tuple(x) for x in g.all_minimal_st_separators())
        min_size_seps = [tuple(x) for x in g.minimum_size_separators()]
        self.failUnless(set(min_size_seps).issubset(min_st_seps))
        self.failUnless(len(set(min_size_seps)) == len(min_size_seps))

        size = len(min_size_seps[0])
        self.failUnless(len(s) != size for s in min_size_seps)
        self.failUnless(sum(1 for s in min_st_seps if len(s) == size) ==
                        len(min_size_seps))


def suite():
    is_separator_suite = unittest.makeSuite(IsSeparatorTests)
    return unittest.TestSuite([is_separator_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())
    
if __name__ == "__main__":
    test()

