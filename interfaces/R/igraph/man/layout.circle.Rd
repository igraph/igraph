\name{layout.circle}
\alias{layout.circle}
\concept{Graph layout}
\title{Graph layout with vertices on a circle.}
\description{
  Place vertices on a circle, in the order of their vertex ids.
}
\usage{
layout.circle (graph, order = V(graph))
}
\arguments{
  \item{graph}{The input graph.}
  \item{order}{The vertices to place on the circle, in the order of
    their desired placement. Vertices that are not included here will
    be placed at (0,0).}
}
\details{
  If you want to order the vertices differently, then permute them using
  the \code{\link{permute.vertices}} function.
}
\value{
  A numeric matrix with two columns, and one row for each vertex.
}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
% \seealso{}
\examples{
## Place vertices on a circle, order them according to their
## community
library(igraphdata)
data(karate)
karate_groups <- optimal.community(karate)
coords <- layout.circle(karate, order =
          order(membership(karate_groups)))
V(karate)$label <- sub("Actor ", "", V(karate)$name)
V(karate)$label.color <- membership(karate_groups)
V(karate)$shape <- "none"
plot(karate, layout = coords)
}
\keyword{graphs}
