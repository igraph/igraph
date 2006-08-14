\name{graph-isomorphism}
\alias{graph.isoclass}
\alias{graph.isocreate}
\alias{graph.isomorphic}
\title{Graph Isomorphism}
\description{These functions deal with graph isomorphism.}
\usage{
graph.isoclass(graph)
graph.isomorphic(graph1, graph2)
graph.isocreate(size, number, directed=TRUE)
}
\arguments{
  \item{graph}{A graph object.}
  \item{graph1,graph2}{Graph objects}
  \item{size}{A numeric integer giving the number of vertices in the
    graph to create. Only three or four are suppported right now.}
  \item{number}{The number of the isomorphism class of the graph to be
    created.}
  \item{directed}{Whether to create a directed graph.}
}
\details{
  \code{graph.isoclass} returns the isomorphism class of a graph, a
  non-negative integer number. Graphs (with the same number of vertices)
  having the same isomorphism class are isomorphic and isomorphic graphs
  always have the same isomorphism class.

  \code{graph.isomorphic} decides whether two graphs are isomorphic.

  \code{graph.isocreate} create a graph from the given isomorphic
  class.
  
  These functions are considered as experimental, as they can
  handle only graphs with three of four vertices.
}
\value{
  \code{graph.isoclass} returns a non-negative integer number.

  \code{graph.isomorphic} returns a logical constant.

  \code{graph.isocreate} returns a graph object.
}
% \references{}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{graph.motifs}}}
\examples{
# create some non-isomorphic graphs
g1 <- graph.isocreate(3, 10)
g2 <- graph.isocreate(3, 11)
graph.isoclass(g1)
graph.isomorphic(g1, g2)
}
\keyword{graphs}
