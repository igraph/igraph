\name{conversion}
\alias{get.adjacency}
\alias{as_adj}
\alias{as_adjacency_matrix}
\alias{get.edgelist}
\alias{as_edgelist}
\concept{Edge list}
\concept{Adjacency list}
\title{Convert a graph to an adjacency matrix or an edge list}
\description{Sometimes it is useful to have a standard representation of
  a graph, like an adjacency matrix or an edge list.}
\usage{
as_adj(graph, type=c("both", "upper", "lower"),
       attr=NULL, edges=FALSE, names=TRUE,
       sparse=getIgraphOpt("sparsematrices"))
as_edgelist(graph, names=TRUE)
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
    (for the multiple edges) is included. This argument is ignored if
    \code{edges} is \code{TRUE}.

    Note that this works only for certain attribute types. If the
    \code{sparse} argumen is \code{TRUE}, then the attribute must be
    either logical or numeric. If the \code{sparse} argument is
    \code{FALSE}, then character is also allowed. The reason for the
    difference is that the \code{Matrix} package does not support
    character sparse matrices yet.
  }
  \item{edges}{Logical scalar, whether to return the edge ids in the
    matrix. For non-existant edges zero is returned.}
  \item{names}{Logical constant.

    For \code{graph.adjacenct} it gives whether to assign row and column names
    to the matrix. These are only assigned if the \code{name} vertex
    attribute is present in the graph.

    for \code{as_edgelist} it gives whether to return a character
    matrix containing vertex names (ie. the \code{name} vertex
    attribute) if they exist or numeric vertex ids.
  }
  \item{sparse}{Logical scalar, whether to create a sparse matrix. The
    \sQuote{\code{Matrix}} package must be installed for creating sparse
    matrices.}
}
\details{
  \code{as_adj} returns the adjacency matrix of a graph, a
  regular \R matrix if \code{sparse} is \code{FALSE}, or a sparse
  matrix, as defined in the \sQuote{\code{Matrix}} package, if
  \code{sparse} if \code{TRUE}.
  
  \code{as_edgelist} returns the list of edges in a graph.  
}
\value{
  A \code{vcount(graph)} by \code{vcount(graph)} (usually) numeric
  matrix for \code{as_adj}. (This can be huge!) Note that a
  non-numeric matrix might be returned if \code{attr} is a non-numeric
  edge attribute.

  A \code{ecount(graph)} by 2 numeric matrix for \code{as_edgelist}.
}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{graph_from_adjacency_matrix}}, \code{\link{read_graph}}}
\examples{
g <- sample_gnp(10, 2/10)
as_edgelist(g)
as_adj(g)
V(g)$name <- letters[1:vcount(g)]
as_adj(g)
E(g)$weight <- runif(ecount(g))
as_adj(g, attr="weight")
}
\keyword{graphs}
