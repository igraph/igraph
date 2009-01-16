\name{graph-isomorphism}
\alias{graph.isoclass}
\alias{graph.isocreate}
\alias{graph.isomorphic}
\alias{graph.isomorphic.vf2}
\alias{graph.count.isomorphisms.vf2}
\alias{graph.count.subisomorphisms.vf2}
\alias{graph.get.isomorphisms.vf2}
\alias{graph.get.subisomorphisms.vf2}
\alias{graph.isoclass.subgraph}
\alias{graph.isomorphic.34}
\alias{graph.isomorphic.bliss}
\alias{graph.subisomorphic.vf2}
\concept{Graph isomorphism}
\concept{Subgraph isomorphism}
\concept{VF2 algorithm}
\concept{BLISS algorithm}
\title{Graph Isomorphism}
\description{These functions deal with graph isomorphism.}
\usage{
graph.isomorphic(graph1, graph2)
graph.isomorphic.34(graph1, graph2)
graph.isomorphic.bliss(graph1, graph2, sh1="fm", sh2="fm")
graph.isomorphic.vf2(graph1, graph2)

graph.count.isomorphisms.vf2(graph1, graph2)
graph.get.isomorphisms.vf2(graph1, graph2)

graph.subisomorphic.vf2(graph1, graph2)
graph.count.subisomorphisms.vf2(graph1, graph2)
graph.get.subisomorphisms.vf2(graph1, graph2)

graph.isoclass(graph)
graph.isoclass.subgraph(graph, vids)
graph.isocreate(size, number, directed=TRUE)
}
\arguments{
  \item{graph}{A graph object.}
  \item{graph1,graph2}{Graph objects}
  \item{size}{A numeric integer giving the number of vertices in the
    graph to create. Only three or four are suppported right now.}
  \item{number}{The number of the isomorphism class of the graph to be
    created.}
  \item{directed}{Whether to create a directed graph.}
  \item{sh1}{Character constant, the heuristics to use in the BLISS
    algorithm, for \code{graph1}. See the \code{sh} argument of
    \code{\link{canonical.permutation}} for possible values.}
  \item{sh2}{Character constant, the heuristics to use in the BLISS
    algorithm, for \code{graph2}. See the \code{sh} argument of
    \code{\link{canonical.permutation}} for possible values.}
  \item{vids}{Numeric vector, the vertex ids of vertices to form the
    induced subgraph for determining the isomorphism class.}
}
\details{
  \code{graph.isomorphic} decides whether two graphs are isomorphic.
  The input graphs must be both directed or both undirected.
  This function is a higher level interface to the other graph
  isomorphism decision functions. Currently it does the following:
  \enumerate{
    \item If the two graphs do not agree in the number of vertices and
    the number of edges then \code{FALSE} is returned.
    \item Otherwise if the graphs have 3 or 4 vertices, then
    \code{igraph.isomorphic.34} is called.
    \item Otherwise if the graphs are directed, then
    \code{igraph.isomorphic.vf2} is called.
    \item Otherwise \code{igraph.isomorphic.bliss} is called.
  }

  \code{igraph.isomorphic.34} decides whether two graphs, both of which
  contains only 3 or 4 vertices, are isomorphic. It works based on a
  precalculated and stored table.

  \code{igraph.isomorphic.bliss} uses the BLISS algorithm by Junttila
  and Kaski, and it works for undirected graphs. For both graphs the
  \code{\link{canonical.permutation}} and then the
  \code{\link{permute.vertices}} function is called to transfer them
  into canonical form; finally the canonical forms are compared.

  \code{graph.isomorphic.vf2} decides whethe two graphs are isomorphic,
  it implements the VF2 algorithm, by Cordella, Foggia et al., see
  references.

  \code{graph.count.isomorphisms.vf2} counts the different isomorphic
  mappings between \code{graph1} and \code{graph2}. (To count
  automorphisms you can supply the same graph twice, but it is better to
  call \code{\link{graph.automorphisms}}.) It uses the VF2 algorithm.

  \code{graph.get.isomorphisms.vf2} calculates all isomorphic mappings
  between \code{graph1} and \code{graph2}. It uses the VF2 algorithm.

  \code{graph.subisomorphic.vf2} decides whether \code{graph2} is
  isomorphic to some subgraph of \code{graph1}. It uses the VF2 algorithm.

  \code{graph.count.subisomorphisms.vf2} counts the different isomorphic
  mappings between \code{graph2} and the subgraphs of \code{graph1}. It
  uses the VF2 algorithm.

  \code{graph.get.subisomorphisms.vf2} calculates all isomorphic
  mappings between \code{graph2} and the subgraphs of \code{graph1}. It
  uses the VF2 algorithm.
  
  \code{graph.isoclass} returns the isomorphism class of a graph, a
  non-negative integer number. Graphs (with the same number of vertices)
  having the same isomorphism class are isomorphic and isomorphic graphs
  always have the same isomorphism class. Currently it can handle only
  graphs with 3 or 4 vertices.

  \code{graph.isoclass.subgraph} calculates the isomorphism class of a
  subgraph of the input graph. Currently it only works for subgraphs
  with 3 or 4 vertices.
  
  \code{graph.isocreate} create a graph from the given isomorphic
  class. Currently it can handle only graphs with 3 or 4 vertices.
}
\note{
  Functions \code{graph.isoclass}, \code{graph.isoclass.subgraph} and
  \code{graph.isocreate} are considered experimental and might be
  reorganized/rewritten later.
}
\value{
  \code{graph.isomorphic} and \code{graph.isomorphic.34} return a
  logical scalar, \code{TRUE} if the input graphs are isomorphic,
  \code{FALSE} otherwise.

  \code{graph.isomorphic.bliss} returns a named list with elements:
  \item{iso}{A logical scalar, whether the two graphs are isomorphic.}
  \item{map12}{A numeric vector, an mapping from \code{graph1} to
    \code{graph2} if \code{iso} is \code{TRUE}, an empty numeric
    vector otherwise.}
  \item{map21}{A numeric vector, an mapping from \code{graph2} to
    \code{graph1} if \code{iso} is \code{TRUE}, an empty numeric
    vector otherwise.}
  \item{info1}{Some information about the canonical form calculation
    for \code{graph1}. A named list, see the return value of
    \code{\link{canonical.permutation}} for details. }
  \item{info2}{Some information about the canonical form calculation
    for \code{graph2}. A named list, see the return value of
    \code{\link{canonical.permutation}} for details. }

  \code{graph.isomorphic.vf2} returns a names list with three elements:
  \item{iso}{A logical scalar, whether the two graphs are isomorphic.}
  \item{map12}{A numeric vector, an mapping from \code{graph1} to
    \code{graph2} if \code{iso} is \code{TRUE}, an empty numeric
    vector otherwise.}
  \item{map21}{A numeric vector, an mapping from \code{graph2} to
    \code{graph1} if \code{iso} is \code{TRUE}, an empty numeric
    vector otherwise.}

  \code{graph.count.isomorphisms.vf2} returns a numeric scalar, an
  integer, the number of isomorphic mappings between the two input
  graphs.

  \code{graph.get.isomorphisms.vf2} returns a list of numeric
  vectors. Every numeric vector is a permutation which takes
  \code{graph2} into \code{graph1}.

  \code{graph.subisomorphic.vf2} returns a named list with three
  elements:
  \item{iso}{Logical scalar, \code{TRUE} if a subgraph of
    \code{graph1} is isomorphic to \code{graph2}.}
  \item{map12}{Numeric vector, empty if \code{iso} is
    \code{FALSE}. Otherwise a mapping from a subgraph of \code{graph1}
    to \code{graph2}. -1 denotes the vertices which are not part of
    the mapping.}
  \item{map21}{Numeric vector, empty if \code{iso} is
    \code{FALSE}. Otherwise a mapping from \code{graph2} into
    \code{graph1}.}

  \code{graph.count.subisomorphisms.vf2} returns a numeric scalar, an
  integer.

  \code{graph.get.subisomorphisms.vf2} returns a list of numeric
  vectors, each numeric vector is an isomorphic mapping from
  \code{graph2} to a subgraph of \code{graph1}.
  
  \code{graph.isoclass} and \code{graph.isoclass.subgraph} return a
  non-negative integer number. 

  \code{graph.isocreate} returns a graph object.
}
\references{
  Tommi Junttila and Petteri Kaski: Engineering an Efficient Canonical
  Labeling Tool for Large and Sparse Graphs, \emph{Proceedings of the
    Ninth Workshop on Algorithm Engineering and Experiments and the
    Fourth Workshop on Analytic Algorithms and Combinatorics.} 2007.

  LP Cordella,  P Foggia, C Sansone, and M Vento:
  An improved algorithm for matching large graphs,
  \emph{Proc. of the 3rd IAPR TC-15 Workshop on Graphbased
    Representations in Pattern Recognition}, 149--159, 2001.
}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{graph.motifs}}}
\examples{
# create some non-isomorphic graphs
g1 <- graph.isocreate(3, 10)
g2 <- graph.isocreate(3, 11)
graph.isoclass(g1)
graph.isoclass(g2)
graph.isomorphic(g1, g2)

# create two isomorphic graphs, by
# permuting the vertices of the first 
g1 <- simplify(barabasi.game(30, m=2, directed=FALSE))
g2 <- permute.vertices(g1, sample(vcount(g1))-1)
# should be TRUE
graph.isomorphic(g1, g2)
graph.isomorphic.bliss(g1, g2)
graph.isomorphic.vf2(g1, g2)
}
\keyword{graphs}
