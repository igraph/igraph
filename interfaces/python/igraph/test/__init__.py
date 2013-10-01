import unittest
from igraph.test import basic, layouts, games, foreign, structural, flow, \
    spectral, attributes, cliques, decomposition, operators, generators, \
    isomorphism, colortests, vertexseq, edgeseq, iterators, bipartite, \
    conversion, rng, separators, indexing, atlas, matching, homepage

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
         rng.suite(),
         separators.suite(),
         indexing.suite(),
         atlas.suite(),
         matching.suite(),
         homepage.suite()] \
    )
    
def run_tests(verbosity=1):
    try:
        # Support for testoob to have nice colored output
        import testoob
        runner = testoob.main
    except ImportError:
        runner = unittest.TextTestRunner(verbosity=verbosity).run
    runner(suite())

# Make nosetest skip run_tests
run_tests.__test__ = False

if __name__ == "__main__":
    run_tests(255)

