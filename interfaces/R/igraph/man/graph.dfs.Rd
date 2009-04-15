\name{graph.dfs}
\alias{graph.dfs}
\title{Depth-first search}
\description{Depth-first search is an algorithm to traverse a graph. It
  starts from a root vertex and tries to go quickly as far from as possible.}
\usage{
graph.dfs (graph, root, neimode = c("out", "in", "all", "total"), 
    order = TRUE, order.out = FALSE, father = FALSE, dist = FALSE, 
    callback = NULL, extra = NULL, rho = parent.frame()) 
}
\arguments{
  \item{graph}{The input graph.}
  \item{root}{The single root vertex to start the search from. }
  \item{neimode}{For directed graphs specifies the type of edges to
    follow. \sQuote{out} follows outgoing, \sQuote{in} incoming
    edges. \sQuote{all} ignores edge directions
    completely. \sQuote{total} is a synonym for \sQuote{all}. This
    argument is ignored for undirected graphs.}
  \item{order}{Logical scalar, whether to return the DFS ordering of the
    vertices.}
  \item{order.out}{Logical scalar, whether to return the ordering based
    on leaving the subtree of the vertex.}
  \item{father}{Logical scalar, whether to return the father of the
    vertices.}
  \item{dist}{Logical scalar, whether to return the distance from the
    root of the search tree.}
  \item{callback}{If not \code{NULL}, then it must be callback
    function. This is called whenever a vertex is visited. See details
    below.}
  \item{extra}{Additional argument to supply to the callback function.}
  \item{rho}{The environment in which the callback function is
    evaluated.}
}
\details{
  The callback function must have the following arguments:
  \describe{
    \item{graph}{The input graph is passed to the callback function
      here.}
    \item{data}{A named numeric vector, with the following entries:
      \sQuote{vid}, the vertex that was just visited and \sQuote{dist},
      its distance from the root of the search tree.}
    \item{extra}{The extra argument.}
  }
  See examples below on how to use the callback function.
}
\value{}
% \references{}
\author{Gabor Csardi \email{Gabor.Csardi@unil.ch} }
\seealso{ \code{\link{igraph.bfs}} for breadth-first search.}
\examples{
}
\keyword{graphs}
