\name{Graphs from adjacency lists}
\alias{graph.adjlist}
\alias{graph_from_adj_list}
\title{Create graphs from adjacency lists}
\description{An adjacency list is a list of numeric vectors, containing
  the neighbor vertices for each vertex. This function creates an igraph
  graph object from such a list.
}
\usage{
graph_from_adj_list(adjlist, mode = c("out", "in", "all", "total"),
              duplicate = TRUE)
}
\arguments{
  \item{adjlist}{The adjacency list. It should be consistent, i.e. the
    maximum throughout all vectors in the list must be less than the
    number of vectors (=the number of vertices in the graph). Note that
    the list is expected to be 0-indexed.}
  \item{mode}{Character scalar, it specifies whether the graph to create
    is undirected (\sQuote{all} or \sQuote{total}) or directed; and in
    the latter case, whether it contains the outgoing (\sQuote{out}) or
    the incoming (\sQuote{in}) neighbors of the vertices.}
  \item{duplicate}{Logical scalar. For undirected graphs it gives
    whether edges are included in the list twice. E.g. if it is
    \code{TRUE} then for an undirected \code{{A,B}} edge
    \code{graph_from_adj_list} expects \code{A} included in the neighbors of
    \code{B} and \code{B} to be included in the neighbors of \code{A}.

    This argument is ignored if \code{mode} is \code{out} or \code{in}.
  }
}
\details{
  Adjacency lists are handy if you intend to do many (small)
  modifications to a graph. In this case adjacency lists are more
  efficient than igraph graphs.

  The idea is that you convert your graph to an adjacency list by
  \code{\link{as_adj_list}}, do your modifications to the graphs and
  finally create again an igraph graph by calling \code{graph_from_adj_list}.
}
\value{
  An igraph graph object.
}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{ \code{\link{as_edgelist}} }
\examples{
## Directed
g <- make_ring(10, dir=TRUE)
al <- as_adj_list(g, mode="out")
g2 <- graph_from_adj_list(al)
graph.isomorphic(g, g2)

## Undirected
g <- make_ring(10)
al <- as_adj_list(g)
g2 <- graph_from_adj_list(al, mode="all")
graph.isomorphic(g, g2)
ecount(g2)
g3 <- graph_from_adj_list(al, mode="all", duplicate=FALSE)
ecount(g3)
which_multiple(g3)
}
\keyword{graphs}
