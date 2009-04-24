\name{graph.bfs}
\alias{graph.bfs}
\title{Breadth-first search}
\description{Breadth-first search is an algorithm to traverse a
  graph. We start from a root vertex and spread along every edge
  \dQuote{simultaneously}. }
\usage{
graph.bfs (graph, root, neimode = c("out", "in", "all", "total"), 
    order = TRUE, rank = FALSE, father = FALSE, pred = FALSE, 
    succ = FALSE, dist = FALSE, callback = NULL, extra = NULL, 
    rho = parent.frame()) 
}
\arguments{
  \item{graph}{The input graph.}
  \item{root}{The single root vertex to start the search from. }
  \item{neimode}{For directed graphs specifies the type of edges to
    follow. \sQuote{out} follows outgoing, \sQuote{in} incoming
    edges. \sQuote{all} ignores edge directions
    completely. \sQuote{total} is a synonym for \sQuote{all}. This
    argument is ignored for undirected graphs.}
  \item{order}{Logical scalar, whether to return the ordering of the
    vertices.}
  \item{rank}{Logical scalar, whether to return the rank of the
    vertices.}
  \item{father}{Logical scalar, whether to return the father of the
    vertices.}
  \item{pred}{Logical scalar, whether to return the predecessors of the
    vertices.}
  \item{succ}{Logical scalar, whether to return the successors of the
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
      \sQuote{vid}, the vertex that was just visited, \sQuote{pred}, its
      predecessor, \sQuote{succ}, its successor, \sQuote{rank}, the rank
      of the current vertex, \sQuote{dist}, its distance from the root
      of the search tree.}
    \item{extra}{The extra argument.}
  }
  See examples below on how to use the callback function.
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
  \item{rank}{Numeric vector. The rank for each vertex.}
  \item{father}{Numeric vector. The father of each vertex, i.e. the
    vertex it was discovered from.}
  \item{pred}{Numeric vector. The previously visited vertex for each
    vertex, or -1 if there was no such vertex.}
  \item{succ}{Numeric vector. The next vertex that was visited after the
    current one, or -1 if there was no such vertex.}
  \item{dist}{Numeric vector, for each vertex its distance from the root
    of the search tree.}
  
  Note that \code{order}, \code{rank}, \code{father}, \code{pred},
  \code{succ} and \code{dist} might be \code{NULL} if their
  corresponding argument is \code{FALSE}, i.e. if their calculation is
  not requested.
}
% \references{}
\author{ Gabor Csardi \email{Gabor.Csardi@unil.ch} }
\seealso{ \code{\link{graph.dfs}} for depth-first search.}
\examples{
## Two rings
graph.bfs(graph.ring(10) \%du\% graph.ring(10), 0, "out",
          TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)

## How to use a callback
f <- function(graph, data, extra) {
  print(data)
  FALSE
}
tmp <- graph.bfs(graph.ring(10) \%du\% graph.ring(10), 0, "out",
                 callback=f)

## How to use a callback to stop the search
## We stop after visiting all vertices in the initial component
f <- function(graph, data, extra) {
 data['succ'] == -1
}
graph.bfs(graph.ring(10) \%du\% graph.ring(10), 0, callback=f)

}
\keyword{graphs}
