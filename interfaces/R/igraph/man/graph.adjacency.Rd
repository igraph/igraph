\name{graph.adjacency}
\alias{graph.adjacency}
\concept{Adjacency matrix}
\concept{Sparse matrix}
\title{Create graphs from adjacency matrices}
\description{\code{graph.adjacency} is a flexible function for creating
  \code{igraph} graphs from adjacency matrices.}
\usage{
graph.adjacency(adjmatrix, mode=c("directed", "undirected", "max",
        "min", "upper", "lower", "plus"), weighted=NULL, diag=TRUE,
        add.colnames=NULL, add.rownames=NA)
}
\arguments{
  \item{adjmatrix}{A square adjacency matrix. From igraph version 0.5.1
    this can be a sparse matrix created with the \code{Matrix} package.}
  \item{mode}{
    Character scalar, specifies how igraph should interpret the supplied
    matrix. See also the \code{weighted} argument, the interpretation
    depends on that too. Possible values are: \code{directed},
    \code{undirected}, \code{upper}, \code{lower}, \code{max},
    \code{min}, \code{plus}. See details below.
  }
  \item{weighted}{This argument specifies whether to create a weighted
    graph from an adjacency matrix. If it is \code{NULL} then an
    unweighted graph is created and the elements of the adjacency matrix
    gives the number of edges between the vertices. If it is a character
    constant then for every non-zero matrix entry an edge is created and
    the value of the entry is added as an edge attribute named by the
    \code{weighted} argument. If it is \code{TRUE} then a weighted graph
    is created and the name of the edge attribute will be
    \code{weight}. See also details below.
  }
  \item{diag}{Logical scalar, whether to include the diagonal of the
    matrix in the calculation. If this is \code{FALSE} then the diagonal
    is zero-d out first.
  }
  \item{add.colnames}{Character scalar, whether to add the column names as
    vertex attributes. If it is \sQuote{\code{NULL}} (the default)
    then, if present, column names are added as vertex attribute
    \sQuote{name}. If \sQuote{\code{NA}} then they will not be added.
    If a character constant, then it gives the name of the vertex
    attribute to add.}
  \item{add.rownames}{Character scalar, whether to add the row names as
    vertex attributes. Possible values the same as the previous
    argument. By default row names are not added. If
    \sQuote{\code{add.rownames}} and \sQuote{\code{add.colnames}} specify the
    same vertex attribute, then the former is ignored. }
}
\details{
  \code{graph.adjacency} creates a graph from an adjacency matrix.

  The order of the vertices are preserved, i.e. the vertex corresponding
  to the first row will be vertex 0 in the graph, etc.

  \code{graph.adjacency} operates in two main modes, depending on the
  \code{weighted} argument.

  If this argument is \code{NULL} then an unweighted graph is
  created and an element of the adjacency matrix gives the number
  of edges to create between the two corresponding vertices.
  The details depend on the value of the \code{mode} argument:
  \itemize{
    \item{\code{directed}}{The graph will be directed and a matrix
      element gives the number of edges between two vertices.}
    \item{\code{undirected}}{This is exactly the same as \code{max},
      for convenience. Note that it is \emph{not} checked whether the
      matrix is symmetric.} 
    \item{\code{max}}{An undirected graph will be created and
      \code{max(A(i,j), A(j,i))} gives the number of edges.}
    \item{\code{upper}}{An undirected graph will be created, only the
      upper right triangle (including the diagonal) is used for the
      number of edges.}
    \item{\code{lower}}{An undirected graph will be created, only the
      lower left triangle (including the diagonal) is used for
      creating the edges.}
    \item{\code{min}}{undirected graph will be created with
      \code{min(A(i,j), A(j,i))} edges between vertex \code{i} and
      \code{j}.}
    \item{\code{plus}}{ undirected graph will be created with
      \code{A(i,j)+A(j,i)} edges between vertex \code{i} and
      \code{j}.}  
  }

  If the \code{weighted} argument is not \code{NULL} then the elements
  of the matrix give the weights of the edges (if they are not zero).
  The details depend on the value of the \code{mode} argument:
  \itemize{
    \item{\code{directed}}{The graph will be directed and a matrix
      element gives the edge weights.}
    \item{\code{undirected}}{First we check that the matrix is
      symmetric. It is an error if not. Then only the upper triangle is
      used to create a weighted undirected graph.}
    \item{\code{max}}{An undirected graph will be created and
      \code{max(A(i,j), A(j,i))} gives the edge weights.}
    \item{\code{upper}}{An undirected graph will be created, only the
      upper right triangle (including the diagonal) is used (for the
      edge weights).}
    \item{\code{lower}}{An undirected graph will be created, only the
      lower left triangle (including the diagonal) is used for
      creating the edges.}
    \item{\code{min}}{An undirected graph will be created,
      \code{min(A(i,j), A(j,i))} gives the edge weights.}
    \item{\code{plus}}{An undirected graph will be created,
      \code{A(i,j)+A(j,i)} gives the edge weights.}
  }
  
}
\value{
  An igraph graph object.
}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\link{graph} and \code{\link{graph.formula}} for
  other ways to create graphs.
}
\examples{
adjm <- matrix(sample(0:1, 100, replace=TRUE, prob=c(0.9,0.1)), ncol=10)
g1 <- graph.adjacency( adjm )
adjm <- matrix(sample(0:5, 100, replace=TRUE,
                      prob=c(0.9,0.02,0.02,0.02,0.02,0.02)), ncol=10)
g2 <- graph.adjacency(adjm, weighted=TRUE)
E(g2)$weight

## various modes for weighted graphs, with some tests
nzs <- function(x) sort(x [x!=0])
adjm <- matrix(runif(100), 10)
adjm[ adjm<0.5 ] <- 0
g3 <- graph.adjacency((adjm + t(adjm))/2, weighted=TRUE,
                      mode="undirected") 

g4 <- graph.adjacency(adjm, weighted=TRUE, mode="max")
all(nzs(pmax(adjm, t(adjm))[upper.tri(adjm)]) == sort(E(g4)$weight))

g5 <- graph.adjacency(adjm, weighted=TRUE, mode="min")
all(nzs(pmin(adjm, t(adjm))[upper.tri(adjm)]) == sort(E(g5)$weight))

g6 <- graph.adjacency(adjm, weighted=TRUE, mode="upper")
all(nzs(adjm[upper.tri(adjm)]) == sort(E(g6)$weight))

g7 <- graph.adjacency(adjm, weighted=TRUE, mode="lower")
all(nzs(adjm[lower.tri(adjm)]) == sort(E(g7)$weight))

g8 <- graph.adjacency(adjm, weighted=TRUE, mode="plus")
d2 <- function(x) { diag(x) <- diag(x)/2; x }
all(nzs((d2(adjm+t(adjm)))[lower.tri(adjm)]) == sort(E(g8)$weight))

g9 <- graph.adjacency(adjm, weighted=TRUE, mode="plus", diag=FALSE)
d0 <- function(x) { diag(x) <- 0 }
all(nzs((d0(adjm+t(adjm)))[lower.tri(adjm)]) == sort(E(g9)$weight))

## row/column names
rownames(adjm) <- sample(letters, nrow(adjm))
colnames(adjm) <- seq(ncol(adjm))
g10 <- graph.adjacency(adjm, weighted=TRUE, add.rownames="code")
summary(g10)
}
\keyword{graphs}
