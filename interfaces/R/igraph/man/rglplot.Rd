\name{rglplot}
\alias{rglplot}
\alias{rglplot.igraph}
\concept{Visualization}
\title{3D plotting of graphs with OpenGL}
\description{Using the \code{rgl} package, \code{rglplot} plots a graph
  in 3D. The plot can be zoomed, rotated, shifted, etc. but the
  coordinates of the vertices is fixed.}
\usage{
rglplot(x, \dots)
}
\arguments{
  \item{x}{The graph to plot.}
  \item{\dots}{Additional arguments, see \code{\link{igraph.plotting}}
    fo the details}
}
\details{
  Note that \code{rglplot} is considered to be highly experimental. It
  is not very useful either. See \code{\link{igraph.plotting}} for the
  possible arguments.
}
\value{\code{NULL}, invisibly.}
% \references{}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{ \code{\link{igraph.plotting}},
  \code{\link{plot.igraph}} for the 2D version,
  \code{\link{tkplot}} for interactive graph drawing in 2D.}
\examples{
\dontrun{
g <- graph.lattice( c(5,5,5) )
coords <- layout.fruchterman.reingold(g, dim=3)
rglplot(g, layout=coords)
}
}
\keyword{graphs}
