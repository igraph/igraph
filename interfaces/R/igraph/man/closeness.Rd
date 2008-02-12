\name{closeness}
\alias{closeness}
\alias{closeness.estimate}
\concept{Closeness centrality}
\title{Closeness centrality of vertices}
\description{Cloness centrality measures how many steps is required
  to access every other vertex from a given vertex.}
\usage{
closeness(graph, v=V(graph), mode = c("all", "out", "in"))
closeness.estimate(graph, vids=V(graph), mode = c("out", "in", "all",
       "total"), cutoff) 
}
\arguments{
  \item{graph}{The graph to analyze.}
  \item{v,vids}{The vertices for which closeness will be calculated.}
  \item{mode}{Character string, defined the types of the paths used for
    measuring the distance in directed graphs. \dQuote{in} measures the
    paths \emph{to} a vertex, \dQuote{out} measures paths \emph{from} a
    vertex, \emph{all} uses undirected paths. This argument is ignored
    for undirected graphs.}
  \item{cutoff}{The maximum path length to consider when calculating the
    betweenness. If zero or negative then there is no such limit.}
}
\details{The closeness centrality of a vertex is defined by the inverse
  of the average length of the shortest paths to/from all the other
  vertices in the graph:

  \deqn{\frac{|V|-1}{\sum_{i\ne v} d_vi}}{(|V|-1)/sum( d(v,i), i != v)}

  If there is no (directed) path between vertex \eqn{v}{\code{v}} and
  \eqn{i}{\code{i}} then the total number of vertices is used in the
  formula instead of the path length.

  \code{closeness.estimate} only considers paths of length
  \code{cutoff} or smaller, this can be run for larger graphs, as the
  running time is not quadratic (if \code{cutoff} is small). If
  \code{cutoff} is zero or negative then the function calculates the
  exact closeness scores.  
}
\value{Numeric vector with the closeness values of all the vertices in
  \code{v}.} 
\references{ Freeman, L.C. (1979). Centrality in Social Networks I:
  Conceptual Clarification. \emph{Social Networks}, 1, 215-239. }
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{betweenness}}, \code{\link{degree}}}
\examples{
g <- graph.ring(10)
g2 <- graph.star(10)
closeness(g)
closeness(g2, mode="in")
closeness(g2, mode="out")
closeness(g2, mode="all")
}
\keyword{graphs}
