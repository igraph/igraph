\name{conversion}
\alias{get.adjacency}
\alias{get.edgelist}
\concept{Edge list}
\concept{Adjacency list}
\title{Convert a graph to an adjacency matrix or an edge list}
\description{Sometimes it is useful to have a standard representation of
  a graph, like an adjacency matrix or an edge list.}
\usage{
get.adjacency(graph, type=c("both", "upper", "lower"),
       attr=NULL, names=TRUE, binary=FALSE, sparse=FALSE)
get.edgelist(graph, names=TRUE)
}
\arguments{
  \item{graph}{The graph to convert.}
  \item{type}{Gives how to create the adjacency matrix for undirected
    graphs. It is ignored for directed graphs. Possible values:
    \code{upper}: the upper right triangle of the matrix is used,
    \code{lower}: the lower left triangle of the matrix is
    used. \code{both}: the whole matrix is used, a symmetric matrix is
    returned.}
  \item{attr}{Either \code{NULL} or a character string giving an edge
    attribute name. If \code{NULL} a traditional adjacency matrix is
    returned. If not \code{NULL} then the values of the given edge
    attribute are included in the adjacency matrix. If the graph has
    multiple edges, the edge attribute of an arbitrarily chosen edge
    (for the multiple edges) is included.}
  \item{names}{Logical constant.

    For \code{graph.adjacenct} it gives whether to assign row and column names
    to the matrix. These are only assigned if the \code{name} vertex
    attribute is present in the graph.

    for \code{get.edgelist} it gives whether to return a character
    matrix containing vertex names (ie. the \code{name} vertex
    attribute) if they exist or numeric vertex ids.
  }
  \item{binary}{Logical, whether to return a binary matrix. This argument
    is ignored if \code{attr} is not \code{NULL}.}
  \item{sparse}{Logical scalar, whether to create a sparse matrix. The
    \sQuote{\code{Matrix}} package must be installed for creating sparse
    matrices.}
}
\details{
  \code{get.adjacency} returns the adjacency matrix of a graph, a
  regular \R matrix if \code{sparse} is \code{FALSE}, or a sparse
  matrix, as defined in the \sQuote{\code{Matrix}} package, if
  \code{sparse} if \code{TRUE}.
  
  \code{get.edgelist} returns the list of edges in a graph.  
}
\value{
  A \code{vcount(graph)} by \code{vcount(graph)} (usually) numeric
  matrix for \code{get.adjacency}. (This can be huge!) Note that a
  non-numeric matrix might be returned if \code{attr} is a non-numeric
  edge attribute.

  A \code{ecount(graph)} by 2 numeric matrix for \code{get.edgelist}.
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{graph.adjacency}}, \code{\link{read.graph}}}
\examples{
g <- erdos.renyi.game(10, 2/10)
get.edgelist(g)
get.adjacency(g)
V(g)$name <- letters[1:vcount(g)]
get.adjacency(g)
E(g)$weight <- runif(ecount(g))
get.adjacency(g, attr="weight")
}
\keyword{graphs}
