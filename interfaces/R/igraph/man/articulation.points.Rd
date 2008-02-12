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
  removal increases the number of connected components in a graph.  
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
g <- add.edges(g, c(which(clu==0), which(clu==1))-1)
articulation.points(g)
}
\keyword{graphs}
