\name{layout}
\alias{layout}
\alias{layout.lgl}
\alias{layout.norm}
\concept{Graph layout}
\title{Generate coordinates for plotting graphs}
\description{Some simple and not so simple functions determining the
  placement of the vertices for drawing a graph.}
\usage{
layout.lgl(graph, \dots, params)
layout.norm(layout, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL,
      zmin = NULL, zmax = NULL)
}
\arguments{
  \item{graph}{The graph to place.}
  \item{params}{The list of function dependent parameters.}
  \item{\dots}{Function dependent parameters, this is an alternative
    notation to the \code{params} argument. For \code{layout.auto} these
    extra parameters are simply passed to the real layout function, if
    one is called.}
  \item{layout}{A matrix with two or three columns, the layout to
    normalize.}
  \item{xmin,xmax}{The limits for the first coordinate, if one of them
    or both are \code{NULL} then no normalization is performed along
    this direction.}
  \item{ymin,ymax}{The limits for the second coordinate, if one of them
    or both are \code{NULL} then no normalization is performed along
    this direction.}
  \item{zmin,zmax}{The limits for the third coordinate, if one of them
    or both are \code{NULL} then no normalization is performed along
    this direction.}  
}
\details{
  These functions calculate the coordinates of the vertices for a graph
  usually based on some optimality criterion.

  \code{layout.lgl} is for large connected graphs, it is similar to the
  layout generator of the Large Graph Layout software
  (\url{http://lgl.sourceforge.net/}). Parameters and
  default values:
  \describe{
    \item{maxiter}{The maximum number of iterations to perform (150).}
    \item{maxdelta}{The maximum change for a vertex during an iteration
      (the number of vertices).}
    \item{area}{The area of the surface on which the vertices are placed
      (square of the number of vertices). }
    \item{coolexp}{The cooling exponent of the simulated annealing (1.5).}
    \item{repulserad}{Cancellation radius for the repulsion (the
      \code{area} times the number of vertices).}
    \item{cellsize}{The size of the cells for the grid. When calculating
      the repulsion forces between vertices only vertices in the same or
      neighboring grid cells are taken into account (the fourth root of
      the number of \code{area}.}
    \item{root}{The id of the vertex to place at the middle of the
      layout. The default value is -1 which means that a random vertex
      is selected.}
  }
  
  \code{layout.norm} normalizes a layout, it linearly transforms each
  coordinate separately to fit into the given limits.

  \code{layout.drl} is another force-driven layout generator, it is
  suitable for quite large graphs. See \code{\link{layout.drl}} for
  details.
}
\value{All these functions return a numeric matrix with at least two
  columns and the same number of lines as the number of vertices.}
%\references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{layout.drl}}, \code{\link{plot.igraph}}, \code{\link{tkplot}}}
% \examples{}
\keyword{graphs}
