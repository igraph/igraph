\name{get.incidence}
\alias{get.incidence}
\concept{Bipartite graph}
\concept{Two-mode network}
\concept{Incidence matrix}
\title{Incidence matrix of a bipartite graph}
\description{This function can return a sparse of dense incidence matrix
of a bipartite network. The incidence matrix is an n times m matrix, n
and m are the number of vertices of the two kinds.x}
\usage{
get.incidence(graph, types=NULL, attr=NULL, names=TRUE, sparse=FALSE)
}
\arguments{
  \item{graph}{The input graph. The direction of the edges is ignored in
  directed graphs.}
  \item{types}{An optional vertex type vector to use instead of the
    \sQuote{\code{type}} vertex attribute. You must supply this argument
    if the graph has no \sQuote{\code{type}} vertex attribute.}
  \item{attr}{Either \code{NULL} or a character string giving an edge
    attribute name. If \code{NULL} a traditional incidence matrix is
    returned. If not \code{NULL} then the values of the given edge
    attribute are included in the incidence matrix. If the graph has
    multiple edges, the edge attribute of an arbitrarily chosen edge
    (for the multiple edges) is included.}
  \item{names}{Logical scalar, if \code{TRUE} and the vertices in the
    graph are named (i.e. the graph has a vertex attribute called
    \sQuote{\code{name}}), then vertex names will be added to the result
    as row and column names. Otherwise the ids of the vertices are used
    as row and column names.}
  \item{sparse}{Logical scalar, if it is \code{TRUE} then a sparse
    matrix is created, you will need the \code{Matrix} package for this.}
}
\details{
  Bipartite graphs have a \sQuote{\code{type}} vertex attribute in
  igraph, this is boolean and \code{FALSE} for the vertices of the first
  kind and \code{TRUE} for vertices of the second kind.

  The vertex ids corresponding to rows and columns in the incidence
  matrix are returned as row/column names.
}
\value{
  A sparse or dense matrix.
}
% \references{}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu} }
\seealso{\code{\link{graph.incidence}} for the opposite operation.}
\examples{
g <- graph.bipartite( c(0,1,0,1,0,0), c(0,1,1,2,2,3) )
get.incidence(g)
}
\keyword{graphs}
