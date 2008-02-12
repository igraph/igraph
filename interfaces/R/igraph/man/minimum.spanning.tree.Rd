\name{minimum.spanning.tree}
\alias{minimum.spanning.tree}
\concept{Minimum spanning tree}
\title{Minimum spanning tree}
\description{A subgraph of a connected graph is a \emph{minimum spanning tree}
  if it is tree, and the sum of its edge weights are the minimal among
  all tree subgraphs of the graph. A minimum spanning forest of a graph
  is the graph consisting of the minimum spanning trees of its
  components.}
\usage{
minimum.spanning.tree(graph, weights=NULL, algorithm=NULL, \dots)
}
\arguments{
  \item{graph}{The graph object to analyze.}
  \item{weights}{Numeric algorithm giving the weights of the edges in
    the graph. The order is determined by the edge ids. This is ignored
    if the \code{unweighted} algorithm is chosen
  }
  \item{algorithm}{The algorithm to use for
    calculation. \code{unweighted} can be used for unwieghted graphs,
    and \code{prim} runs Prim's algorithm for weighted graphs.
    If this is \code{NULL} then igraph tries to select the algorithm
    automatically: if the graph has an edge attribute called
    \code{weight} of the \code{weights} argument is not \code{NULL} then
    Prim's algorithm is chosen, otherwise the unwweighted algorithm is
    performed.
  }
  \item{\dots}{Additional arguments, unused.}
}
\details{
  If the graph is unconnected a minimum spanning forest is returned.
}
\value{
  A graph object with the minimum spanning forest. (To check that it is
  a tree check that the number of its edges is \code{vcount(graph)-1}.)
  The edge and vertex attributes of the original graph are preserved in
  the result.
}
\references{ Prim, R.C. 1957. Shortest connection networks and some
  generalizations \emph{Bell System Technical Journal}, 37 1389--1401.}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{clusters}}}
\examples{
g <- erdos.renyi.game(100, 3/100)
mst <- minimum.spanning.tree(g)
}
\keyword{graphs}
