\name{minimum.spanning.tree}
\alias{minimum.spanning.tree}
\title{Minimum spanning tree}
\description{A subgraph of a connected graph is a minimum spanning tree
  if it is tree, and the sum of its edge weights are the minimal among
  all tree subgraphs of the graph. A minimum spanning forest of a graph
  is the graph consisting of the minimum spanning trees of its
  components.}
\usage{
minimum.spanning.tree(graph, weights=NULL, algorithm="prim", \dots)
}
\arguments{
  \item{graph}{The graph object to analyze.}
  \item{weights}{Logical. If \code{TRUE}, the spanning tree is
    calculated based of the weights given in the \code{weight} edge
    attribute. If \code{FALSE} every edge is treated as having weight
    one. By default (\code{NULL}) this option is set to \code{TRUE} if
    the graph has an edge attribute called \code{weight} and to
    \code{FALSE} otherwise.
  }
  \item{algorithm}{The algorithm to use for the weighted
    calculation. For unweighted graphs a simple breath-first search is
    used, always. Right now only Prim's algorithm is implemented
    (\code{prim}), but this gives satisfactory results.}
  \item{\dots}{Additional arguments to be passed to the graph
    constructor.}
}
\details{
  Some remarks.

  The calculation is quite slow for adjacency lists, because of the poor
  implementation of edge attributes.

  If the graph is unconnected a minimum spanning forest is returned.
}
\value{
  A graph object with the minimum spanning forest. (To check that it is
  a tree check that the number of its edges is \code{vcount(graph)-1}.
}
\references{ Prim, R.C. 1957. ``Shortest connection networks and some
  generalizations'' Bell System Technical Journal, 37 1389--1401.}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{clusters}}}
\examples{
g <- erdos.renyi.game(100, 3/100)
mst <- minimum.spanning.tree(g)
}
\keyword{graphs}
