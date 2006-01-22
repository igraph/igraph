\name{conversion}
\alias{get.adjacency}
\alias{get.edgelist}
\title{Convert a graph to an adjacency matrix or an edge list}
\description{Sometimes it is useful to have a standard representation of
  a graph, like an adjacency matrix or an edge list.}
\usage{
get.adjacency(graph, type="both")
get.edgelist(graph)
}
\arguments{
  \item{graph}{The graph to convert.}
  \item{type}{Gives how to create the adjacency matrix for undirected
    graphs. It is ignored for directed graphs. Possible values:
    \code{upper}: the upper right triangle of the matrix is used,
    \code{lower}: the lower left triangle of the matrix is
    used. \code{both}: the whole matrix is used, a symmetric matrix is
    returned.} 
}
\details{
  \code{get.adjacency} returns the adjacency matrix of a graph.

  \code{get.edgelist} returns the list of edges in a graph.  
}
\value{
  A \code{vcount(graph)} by \code{vcount(graph)} numeric matrix for
  \code{get.adjacency}. (This can be huge!)

  A \code{ecount(graph)} by 2 numeric matrix for \code{get.edgelist}.
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{graph.adjacency}}, \code{\link{read.graph}}}
\examples{
g <- erdos.renyi.game(10, 2/10)
get.adjacency(g)
get.edgelist(g)
}
\keyword{graphs}
