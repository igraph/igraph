import unittest
from igraph.test import basic, layouts, games, foreign, structural, flow, \
    spectral, attributes, cliques, decomposition, operators, generators, \
    isomorphism, colortests, vertexseq, edgeseq, iterators, bipartite, \
    conversion, rng

def suite():
    return unittest.TestSuite( \
        [basic.suite(),
         layouts.suite(),
         generators.suite(),
         games.suite(),
         foreign.suite(),
         structural.suite(),
         flow.suite(),
         spectral.suite(),
         attributes.suite(),
         vertexseq.suite(),
         edgeseq.suite(),
         cliques.suite(),
         decomposition.suite(),
         conversion.suite(),
         operators.suite(),
         isomorphism.suite(),
         iterators.suite(),
         bipartite.suite(),
         colortests.suite(),
         rng.suite()] \
    )
    
def test():
    try:
        # Support for testoob to have nice colored output
        import testoob
        testoob.main(suite())
    except ImportError:
        runner = unittest.TextTestRunner(verbosity=2)
        runner.run(suite())

if __name__ == "__main__": test()

