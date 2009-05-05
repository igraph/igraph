\name{diameter}
\alias{diameter}
\alias{get.diameter}
\alias{farthest.nodes}
\concept{Diameter}
\title{Diameter of a graph}
\description{The diameter of a graph is the length of the longest
  geodesic.}
\usage{
diameter(graph, directed = TRUE, unconnected = TRUE, weights = NULL)
get.diameter (graph, directed = TRUE, unconnected = TRUE, weights = NULL) 
farthest.nodes (graph, directed = TRUE, unconnected = TRUE, weights = NULL) 
}
\arguments{
  \item{graph}{The graph to analyze.}
  \item{directed}{Logical, whether directed or undirected paths are to
    be considered. This is ignored for undirected graphs.}
  \item{unconnected}{Logical, what to do if the graph is unconnected. If
    FALSE, the function will return the largest possible diameter which
    is the number of vertices. If TRUE, the diameters of the connected
    components will be calculated and the largest one will be returned.}
  \item{weights}{Optional positive weight vector for calculating
    weighted distances. If the graph has a \code{weight} edge
    attribute, then this is used by default.}
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

## Weighted diameter
set.seed(1)
g <- graph.ring(10)
E(g)$weight <- sample(seq_len(ecount(g)))
diameter(g)
get.diameter(g)
diameter(g, weights=NA)
get.diameter(g, weights=NA)
}
\keyword{graphs}
