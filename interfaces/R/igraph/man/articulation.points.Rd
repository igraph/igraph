\name{articulation_points}
\alias{articulation.points}
\alias{articulation_points}
\concept{Articulation point}
\title{Articulation points of a graph}
\description{Articuation points or cut vertices are vertices whose
  removal increases the number of connected components in a graph.}
\usage{
articulation_points(graph)
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
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{biconnected_components}}, \code{\link{components}},
  \code{\link{is_connected}}, \code{\link{vertex_connectivity}}}
\examples{
g <- disjoint_union( make_full_graph(5), make_full_graph(5) )
clu <- components(g)$membership
g <- add_edges(g, c(match(1, clu), match(2, clu)) )
articulation_points(g)
}
\keyword{graphs}
