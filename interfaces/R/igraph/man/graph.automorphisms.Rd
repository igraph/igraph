\name{graph.automorphisms}
\alias{graph.automorphisms}
\concept{Graph automorphism}
\title{Number of automorphisms}
\description{Calculate the number of automorphisms of a graph,
  i.e. the number of isomorphisms to itself.}
\usage{
graph.automorphisms(graph, sh="fm")
}
\arguments{
  \item{graph}{The input graph, it is treated as undirected.}
  \item{sh}{The splitting heuristics for the BLISS algorithm. Possible
    values are: \sQuote{\code{f}}: first non-singleton cell,
    \sQuote{\code{fl}}: first largest non-singleton cell,
    \sQuote{\code{fs}}: first smallest non-singleton cell,
    \sQuote{\code{fm}}: first maximally non-trivially connected
    non-singleton cell, \sQuote{\code{flm}}: first largest maximally
    non-trivially connected non-singleton cell, \sQuote{\code{fsm}}:
    first smallest maximally non-trivially connected non-singleton
    cell.}
}
\details{
  An automorphism of a graph is a permutation of its vertices which
  brings the graph into itself.

  This function calculates the number of automorphism of a graph using
  the BLISS algorithm. See also the BLISS homepage at 
  \url{http://www.tcs.hut.fi/Software/bliss/index.html}.
}
\value{
  A named list with the following members:
  \item{group_size}{The size of the automorphism group of the
    input graph, as a string. This number is exact if igraph was
    compiled with the GMP library, and approximate otherwise.}
  \item{nof_nodes}{The number of nodes in the search tree.}
  \item{nof_leaf_nodes}{The number of leaf nodes in the search
    tree.}
  \item{nof_bad_nodes}{Number of bad nodes.}
  \item{nof_canupdates}{Number of canrep updates.}
  \item{max_level}{Maximum level.}
}
\references{
  Tommi Junttila and Petteri Kaski: Engineering an Efficient Canonical
  Labeling Tool for Large and Sparse Graphs, \emph{Proceedings of the
    Ninth Workshop on Algorithm Engineering and Experiments and the
    Fourth Workshop on Analytic Algorithms and Combinatorics.} 2007.
}
\author{ Tommi Juntilla \email{google@for.it} for BLISS and Gabor Csardi
  \email{csardi@rmki.kfki.hu} for the igraph glue code and this manual
  page.} 
\seealso{\code{\link{canonical.permutation}}, \code{\link{permute.vertices}}}
\examples{
## A ring has n*2 automorphisms, you can "turn" it by 0-9 vertices
## and each of these graphs can be "flipped"
g <- graph.ring(10)
graph.automorphisms(g)
}
\keyword{graphs}
