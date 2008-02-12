\name{shortest.paths}
\alias{shortest.paths}
\alias{get.shortest.paths}
\alias{get.all.shortest.paths}
\alias{average.path.length}
\alias{path.length.hist}
\concept{Shortest path}
\concept{Geodesic}
\title{Shortest (directed or undirected) paths between vertices}
\description{\code{shortest.paths} calculates the length of all the
  shortest paths from or to the vertices in the
  network. \code{get.shortest.paths} calculates one shortest path (the
  path itself, and not just its length) from or to the given vertex.}
\usage{
shortest.paths(graph, v=V(graph), mode = c("all", "out", "in"))
get.shortest.paths(graph, from, to=V(graph), mode = c("all", "out", "in"))
get.all.shortest.paths(graph, from, to = V(graph), mode = c("all", "out", "in")) 
average.path.length(graph, directed=TRUE, unconnected=TRUE)
path.length.hist (graph, directed = TRUE, verbose = igraph.par("verbose")) 
}
\arguments{
  \item{graph}{The graph to work on.}
  \item{v}{Numeric vector, the vertices from or to which the shortest
    paths will be calculated.}
  \item{mode}{Character constant, gives whether the shortest paths to or
    from the given vertices should be calculated for directed graphs. If
    \code{out} then the shortest paths \emph{from} the vertex, if
    \code{in} then \emph{to} it will be considered. If \code{all}, the
    default, then the corresponding undirected graph will be used,
    ie. not directed paths are searched. This argument is ignored for
    undirected graphs.}
  \item{from}{Numeric constant, the vertex from or to the shortest paths
    will be calculated. Note that right now this is not a vector of
    vertex ids, but only a single vertex.}
  \item{to}{Numeric vector, only the shortest paths to these vertices
    will be calculated. Defaults to all vertices.}
  \item{directed}{Whether to consider directed paths in directed graphs,
    this argument is ignored for undirected graphs.}
  \item{unconnected}{What to do if the graph is unconnected (not
    strongly connected if directed paths are considered). If TRUE only
    the lengths of the existing paths are considered and averaged; if
    FALSE the length of the missing paths are counted having length
    \code{vcount(graph)}, one longer than the longest possible geodesic
    in the network.}
  \item{verbose}{Logical scalar, whether to draw a progress meter while
    the calculation is running.}
}
\details{The shortest paths (also called geodesics) are calculated by
  using breath-first search in the graph. Right now edge weights are not
  used, ie. every edge weight is one.

  Note that \code{shortest.paths} is able to calculate the path length
  from or to many vertices at the same time, but
  \code{get.shortest.paths} works from one source only. This might
  change in the future.

  Also note that \code{get.shortest.paths} gives only one shortest path,
  however more than one might exists between two vertices.

  \code{get.all.shortest.paths} calculates all shortest paths from a
  vertex to other vertices given in the \code{to} argument.

  \code{path.length.hist} calculates a histogram, by calculating the
  shortest path length between each pair of vertices. For directed
  graphs both directions are considered, so every pair of vertices
  appears twice in the histogram.
}
\value{
  For \code{shortest.paths} a numeric matrix with \code{vcount(graph)}
  columns and \code{length(v)} rows. The shortest path length from a vertex to
  itself is always zero.

  For \code{get.shortest.paths} a list of length
  \code{vcount(graph)}. List element \code{i} contains the vertex ids on
  the path from vertex \code{from} to vertex \code{i} (or the other way
  for directed graphs depending on the \code{mode} argument). The vector
  also contains \code{from} and \code{i} as the first and last
  elements. If \code{from} is the same as \code{i} then it is only
  included once. If there is no path between two vertices then a
  numeric vector of length zero is returned as the list element.

  For \code{get.all.shortest.paths} a list is returned, each list
  element contains a shortest path from \code{from} to a vertex in
  \code{to}. The shortest paths to the same vertex are collected into
  consecutive elements of the list.
  
  For \code{average.path.length} a single number is returned.

  \code{path.length.hist} returns a named list with two entries:
  \code{res} is a numeric vector, the histogram of distances,
  \code{unconnected} is a numeric scalar, the number of pairs for which
  the first vertex is not reachable from the second. The sum of the two
  entries is always \eqn{n(n-1)} for directed graphs and \eqn{n(n-1)/2}
  for undirected graphs.
}
\references{ West, D.B. (1996). \emph{Introduction to Graph Theory.} Upper
  Saddle River, N.J.: Prentice Hall. }
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
%\seealso{}
\examples{
g <- graph.ring(10)
shortest.paths(g)
get.shortest.paths(g, 5)
get.all.shortest.paths(g, 0, 5:7)
average.path.length(g)
}
\keyword{graphs}
