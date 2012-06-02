\name{graph.dfs}
\alias{graph.dfs}
\title{Depth-first search}
\description{Depth-first search is an algorithm to traverse a graph. It
  starts from a root vertex and tries to go quickly as far from as possible.}
\usage{
graph.dfs (graph, root, neimode = c("out", "in", "all", "total"), 
    unreachable = TRUE, order = TRUE, order.out = FALSE,
    father = FALSE, dist = FALSE, in.callback = NULL,
    out.callback = NULL, extra = NULL, rho = parent.frame()) 
}
\arguments{
  \item{graph}{The input graph.}
  \item{root}{The single root vertex to start the search from. }
  \item{neimode}{For directed graphs specifies the type of edges to
    follow. \sQuote{out} follows outgoing, \sQuote{in} incoming
    edges. \sQuote{all} ignores edge directions
    completely. \sQuote{total} is a synonym for \sQuote{all}. This
    argument is ignored for undirected graphs.}
  \item{unreachable}{Logical scalar, whether the search should visit
    the vertices that are unreachable from the given root
    vertex (or vertices). If \code{TRUE}, then additional searches are
    performed until all vertices are visited.}
  \item{order}{Logical scalar, whether to return the DFS ordering of the
    vertices.}
  \item{order.out}{Logical scalar, whether to return the ordering based
    on leaving the subtree of the vertex.}
  \item{father}{Logical scalar, whether to return the father of the
    vertices.}
  \item{dist}{Logical scalar, whether to return the distance from the
    root of the search tree.}
  \item{in.callback}{If not \code{NULL}, then it must be callback
    function. This is called whenever a vertex is visited. See details
    below.}
  \item{out.callback}{If not \code{NULL}, then it must be callback
    function. This is called whenever the subtree of a vertex is
    completed by the algorithm. See details below.}
  \item{extra}{Additional argument to supply to the callback function.}
  \item{rho}{The environment in which the callback function is
    evaluated.}
}
\details{
  The callback functions must have the following arguments:
  \describe{
    \item{graph}{The input graph is passed to the callback function
      here.}
    \item{data}{A named numeric vector, with the following entries:
      \sQuote{vid}, the vertex that was just visited and \sQuote{dist},
      its distance from the root of the search tree.}
    \item{extra}{The extra argument.}
  }
  See examples below on how to use the callback functions.
}
\value{
  A named list with the following entries:
  \item{root}{Numeric scalar. The root vertex that was used as the
    starting point of the search.}
  \item{neimode}{Character scalar. The \code{neimode} argument of the
    function call. Note that for undirected graphs this is always
    \sQuote{all}, irrespectively of the supplied value.}
  \item{order}{Numeric vector. The vertex ids, in the order in which
    they were visited by the search.}
  \item{order.out}{Numeric vector, the vertex ids, in the order of the
    completion of their subtree.}
  \item{father}{Numeric vector. The father of each vertex, i.e. the
    vertex it was discovered from.}
  \item{dist}{Numeric vector, for each vertex its distance from the root
    of the search tree.}

  Note that \code{order}, \code{order.out}, \code{father}, and
  \code{dist} might be \code{NULL} if their corresponding argument is
  \code{FALSE}, i.e. if their calculation is not requested.  
}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com} }
\seealso{ \code{\link{graph.bfs}} for breadth-first search.}
\examples{
## A graph with two separate trees
graph.dfs(graph.tree(10) \%du\% graph.tree(10), root=1, "out",
          TRUE, TRUE, TRUE, TRUE)

## How to use a callback
f.in <- function(graph, data, extra) {
  cat("in:", paste(collapse=", ", data), "\n")
  FALSE
}
f.out <- function(graph, data, extra) {
  cat("out:", paste(collapse=", ", data), "\n")
  FALSE
}
tmp <- graph.dfs(graph.tree(10), root=1, "out",
                 in.callback=f.in, out.callback=f.out)

## Terminate after the first component, using a callback
f.out <- function(graph, data, extra) {
 data['vid'] == 1
}
tmp <- graph.dfs(graph.tree(10) \%du\% graph.tree(10), root=1,
                 out.callback=f.out)

}
\keyword{graphs}
