\name{diameter}
\alias{diameter}
\alias{get.diameter}
\alias{farthest.nodes}
\title{Diameter of a graph}
\description{The diameter of a graph is the length of the longest
  geodesic.}
\usage{
diameter(graph, directed = TRUE, unconnected = TRUE)
get.diameter (graph, directed = TRUE, unconnected = TRUE) 
farthest.nodes (graph, directed = TRUE, unconnected = TRUE) 
}
\arguments{
  \item{graph}{The graph to analyize.}
  \item{directed}{Logical, whether directed or undirected paths are to
    be considered. This is ignored for undirected graphs.}
  \item{unconnected}{Logical, what to do if the graph is unconnected. If
    FALSE, the function will return the largest possible diameter which
    is the number of vertices. If TRUE, the diameters of the connected
    components will be calculated and the largest one will be returned.}
}
\details{The diameter is calculated by using a breadth-first search like
  method.

  \code{get.diameter} returns a path with the actual diameter. If there
  are many shortest paths of the length of the diameter, then it returns
  the first one found.

  \code{farthest.points} returns two vertex ids, the vertices which are
  connected by the diameter path.
}
\value{A numeric constant for \code{diameter}, a numeric vector for
  \code{get.diameter} and a numeric vector of length two for
  \code{farthest.nodes}.
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{shortest.paths}}}
\examples{
g <- graph.ring(10)
g2 <- delete.edges(g, c(0,1,0,9))
diameter(g2, unconnected=TRUE)
diameter(g2, unconnected=FALSE)
}
\keyword{graphs}
