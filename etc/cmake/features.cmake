include(helpers)

include(tls)
include(lto)

option(IGRAPH_GLPK_SUPPORT "Compile igraph with GLPK support" ON)
tristate(IGRAPH_GRAPHML_SUPPORT "Compile igraph with GraphML support" AUTO)
tristate(IGRAPH_OPENMP_SUPPORT "Use OpenMP for parallelization" AUTO)
