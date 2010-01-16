\name{layout.star}
\alias{layout.star}
\concept{Graph layout}
\title{Generate coordinates to place the vertices of a graph in a
  star-shape}
\description{
  A simple layout generator, that places one vertex in the center of a
  circle and the rest of the vertices equidistantly on the perimeter.
}
\usage{
layout.star(graph, center = V(graph)[0], order = NULL)
}
\arguments{
  \item{graph}{The graph to layout.}
  \item{center}{The id of the vertex to put in the center. By default it
    is the first vertex.}
  \item{order}{Numeric vector, the order of the vertices along the
    perimeter. The default ordering is given by the vertex ids.}
}
\details{
  It is possible to choose the vertex that will be in the center, and
  the order of the vertices can be also given.
}
\value{A matrix with two columns and as many rows as the number of
  vertices in the input graph.}
%\references{}
\author{Gabor Csardi \email{Gabor.Csardi@unil.ch}}
\seealso{\code{\link{layout}} and \code{\link{layout.drl}} for other
  layout algorithms, \code{\link{plot.igraph}} and \code{\link{tkplot}}
  on how to plot graphs and \code{\link{graph.star}} on how to create
  ring graphs.}
\examples{
g <- graph.star(10)
layout.star(g)
}
\keyword{graphs}
