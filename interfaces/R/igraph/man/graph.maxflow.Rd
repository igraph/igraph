\name{graph.maxflow}
\alias{graph.maxflow}
\alias{graph.mincut}
\concept{Maximum flow}
\concept{Minimum cut}
\title{Maximum flow in a network}
\description{In a graph where each edge has a given flow capacity the
  maximal flow between two vertices is calculated.}
\usage{
graph.maxflow(graph, source, target, capacity=NULL)
graph.mincut(graph, source=NULL, target=NULL, capacity=NULL,
       value.only = TRUE)
}
\arguments{
  \item{graph}{The input graph.}
  \item{source}{The id of the source vertex.}
  \item{target}{The id of the target vertex (sometimes also called sink).}
  \item{capacity}{Vector giving the capacity of the edges. If this is
    \code{NULL} (the default) then the \code{capacity} edge attribute is
    used.}
  \item{value.only}{Logical scalar, if \code{TRUE} only the minumum cut
    value is returned, if \code{FALSE} the edges in the cut and a the
    two (or more) partitions are also returned. This currently only
    works for undirected graphs.
  }
}
\details{
  \code{graph.maxflow} calculates the maximum flow between two vertices
  in a weighted (ie. valued) graph. A flow from \code{source} to
  \code{target} is an assignment of non-negative real numbers to the
  edges of the graph, satisfying two properties: (1) for each edge the
  flow (ie. the assigned number) is not more than the capacity of the
  edge (the \code{capacity} parameter or edge attribute), (2) for every
  vertex, except the source and the target the incoming flow is the same
  as the outgoing flow. The value of the flow is the incoming flow of
  the \code{target} vertex. The maximum flow is the flow of maximum
  value.

  \code{graph.mincut} calculates the minimum st-cut between two vertices
  in a graph (if the \code{source} and \code{target} arguments are
  given) or the minimum cut of the graph (if both \code{source} and
  \code{target} are \code{NULL}).

  The minimum st-cut between \code{source} and \code{target} is the
  minimum total weight of edges needed to remove to eliminate all paths from
  \code{source} to \code{target}.

  The minimum cut of a graph is the minimum total weight of the edges
  needed to remove to separate the graph into (at least) two
  components. (Which is to make the graph \emph{not} strongly connected
  in the directed case.)
  
  The maximum flow between two vertices in a graph is the same as the minimum
  st-cut, so \code{graph.maxflow} and \code{graph.mincut} essentially
  calculate the same quantity, the only difference is that
  \code{graph.mincut} can be invoked without giving the \code{source}
  and \code{target} arguments and then minimum of all possible minimum
  cuts is calculated.

  Starting from version 0.5 igraph can return the edges in the minimum
  cut for undirected graphs.
}
\value{
  A numeric constant, the maximum flow, or the minimum cut, except if
  \code{value.only=FALSE} for \code{graph.mincut}. In this case a named
  list with components:
  \item{value}{Numeric scalar, the cut value.}
  \item{cut}{Numeric vector, the edges in the cut.}
  \item{partition1}{The vertices in the first partition after the cut
    edges are removed. Note that these vertices might be actually in
    different components (after the cut edges are removed), as the
    graph may fall apart into more than two components.}
  \item{partition2}{The vertices in the second partition after the cut
    edges are removed. Note that these vertices might be actually in
    different components (after the cut edges are removed), as the
    graph may fall apart into more than two components.}
}
\references{
  A. V. Goldberg and R. E. Tarjan: ``A New Approach to the Maximum Flow
  Problem'' Journal of the ACM 35:921-940, 1988.  
}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{shortest.paths}}, \code{\link{edge.connectivity}},
  \code{\link{vertex.connectivity}}}
\examples{
g <- graph.ring(100)
graph.mincut(g, capacity=rep(1,vcount(g)))
graph.mincut(g, value.only=FALSE, capacity=rep(1,vcount(g)))
}
\keyword{graphs}
