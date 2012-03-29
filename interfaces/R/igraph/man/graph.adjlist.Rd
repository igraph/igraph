\name{Graphs from adjacency lists}
\alias{graph.adjlist}
\title{Create graphs from adjacency lists}
\description{An adjacency list is a list of numeric vectors, containing
  the neighbor vertices for each vertex. This function creates an igraph
  graph object from such a list.
}
\usage{
graph.adjlist(adjlist, directed = TRUE, duplicate = TRUE)
}
\arguments{
  \item{adjlist}{The adjacency list. It should be consistent, i.e. the
    maximum throughout all vectors in the list must be less than the
    number of vectors (=the number of vertices in the graph). Note that
    the list is expected to be 0-indexed.}
  \item{directed}{Logical scalar, whether or not to create a directed
    graph.}
  \item{duplicate}{Logical scalar. For undirected graphs it gives
    whether edges are included in the list twice. E.g. if it is
    \code{TRUE} then for an undirected \code{{A,B}} edge
    \code{graph.adjlist} expects \code{A} included in the neighbors of
    \code{B} and \code{B} to be included in the neighbors of \code{A}.

    This argument is ignored if \code{directed} is \code{TRUE}.
  }
}
\details{
  Adjacency lists are handy if you intend to do many (small)
  modifications to a graph. In this case adjacency lists are more
  efficient than igraph graphs.

  The idea is that you convert your graph to an adjacency list by
  \code{\link{get.adjlist}}, do your modifications to the graphs and
  finally create again an igraph graph by calling \code{graph.adjlist}.
}
\value{
  An igraph graph object.
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{ \code{\link{get.edgelist}} }
\examples{
## Directed
g <- graph.ring(10, directed=TRUE)
al <- get.adjlist(g, mode="out")
g2 <- graph.adjlist(al)
graph.isomorphic(g, g2)

## Undirected
g <- graph.ring(10)
al <- get.adjlist(g)
g2 <- graph.adjlist(al, directed=FALSE)
graph.isomorphic(g, g2)
ecount(g2)
g3 <- graph.adjlist(al, directed=FALSE, duplicate=FALSE)
ecount(g3)
is.multiple(g3)
}
\keyword{graphs}
