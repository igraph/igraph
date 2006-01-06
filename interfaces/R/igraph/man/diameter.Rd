\name{diameter}
\alias{diameter}
\title{Diameter of a graph}
\description{The diameter of a graph is the length of the longest
  geodesic.}
\usage{
diameter(graph, directed = TRUE, unconnected = TRUE)
}
\arguments{
  \item{graph}{The graph to analyize.}
  \item{directed}{Logical, whether directed or undirected paths are to
    be considered. This is ignored for undirected graphs.}
  \item{unconnected}{Logical, what to do if the graph is unconnected. If
    FALSE, the function will return the largest possible diameter which
    is the number of vertices. If TRUE, the diameter within connected
    clusters will be considered only.}
}
\details{The diameter is calculated by using a breadth-first search like
  method.}
\value{A numeric constant.}
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
