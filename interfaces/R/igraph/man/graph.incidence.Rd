\name{graph.incidence}
\alias{graph.incidence}
\concept{Bipartite graph}
\concept{Two-mode network}
\title{Create graphs from an incidence matrix}
\description{\code{graph.incidence} creates a bipartite igraph graph
  from an incidence matrix.}
\usage{
graph.incidence(incidence, directed = FALSE, mode = c("all", "out", 
    "in", "total"), multiple = FALSE, weighted = NULL, add.names = NULL)
}
\arguments{
  \item{incidence}{The input incidence matrix. It can also be a sparse
    matrix from the \code{Matrix} package.}
  \item{directed}{Logical scalar, whether to create a directed graph.}
  \item{mode}{A character constant, defines the direction of the edges
    in directed graphs, ignored for undirected graphs. If
    \sQuote{\code{out}}, then edges go from vertices of the first kind
    (corresponding to rows in the incidence matrix) to vertices of the
    second kind (columns in the incidence matrix). If
    \sQuote{\code{in}}, then the opposite direction is used. If
    \sQuote{\code{all}} or \sQuote{\code{total}}, then mutual edges are
    created. }
  \item{multiple}{Logical scalar, specifies how to interpret the matrix
    elements. See details below.}
  \item{weighted}{This argument specifies whether to create a weighted
    graph from the incidence matrix. If it is \code{NULL} then an
    unweighted graph is created and the \code{multiple} argument is used
    to determine the edges of the graph. If it is a character
    constant then for every non-zero matrix entry an edge is created and
    the value of the entry is added as an edge attribute named by the
    \code{weighted} argument. If it is \code{TRUE} then a weighted graph
    is created and the name of the edge attribute will be
    \sQuote{\code{weight}}.
  }
  \item{add.names}{A character constant, \code{NA} or
    \code{NULL}. \code{graph.incidence} can add the row and column names
    of the incidence matrix as vertex attributes. If this argument is
    \code{NULL} (the default) and the incidence matrix has both row and
    column names, then these are added as the \sQuote{\code{name}} vertex
    attribute. If you want a different vertex attribute for this, then
    give the name of the attributes as a character string. If this
    argument is \code{NA}, then no vertex attributes (other than type)
    will be added.
  }
}
\details{
  Bipartite graphs have a \sQuote{\code{type}} vertex attribute in
  igraph, this is boolean and \code{FALSE} for the vertices of the first
  kind and \code{TRUE} for vertices of the second kind.

 \code{graph.incidence} can operate in two modes, depending on the
 \code{multiple} argument. If it is \code{FALSE} then a single edge is
 created for every non-zero element in the incidence matrix. If
 \code{multiple} is \code{TRUE}, then the matrix elements are rounded up
 to the closest non-negative integer to get the number of edges to
 create between a pair of vertices.
}
\value{
  A bipartite igraph graph. In other words, an igraph graph that has a
  vertex attribute \code{type}. 
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{graph.bipartite}} for another way to create
  bipartite graphs}
\examples{
inc <- matrix(sample(0:1, 15, repl=TRUE), 3, 5)
colnames(inc) <- letters[1:5]
rownames(inc) <- LETTERS[1:3]
graph.incidence(inc)
}
\keyword{graphs}
