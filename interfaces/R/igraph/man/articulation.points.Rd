\name{articulation.points}
\alias{articulation.points}
\concept{Articulation point}
\title{Articulation points of a graph}
\description{Articuation points or cut vertices are vertices whose
  removal increases the number of connected components in a graph.}
\usage{
articulation.points(graph)
}
\arguments{
  \item{graph}{The input graph. It is treated as an undirected graph,
    even if it is directed.}
}
\details{
  Articuation points or cut vertices are vertices whose
  removal increases the number of connected components in a graph. If
  the original graph was connected, then the removal of a single
  articulation point makes it undirected. If a graph contains no
  articulation points, then its vertex connectivity is at least two.  
}
\value{A numeric vector giving the vertex ids of the articulation points
  of the input graph.}
%\references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{biconnected.components}}, \code{\link{clusters}},
  \code{\link{is.connected}}, \code{\link{vertex.connectivity}}}
\examples{
g <- graph.disjoint.union( graph.full(5), graph.full(5) )
clu <- clusters(g)$membership
g <- add.edges(g, c(match(1, clu), match(2, clu)) )
articulation.points(g)
}
\keyword{graphs}
