\name{layout.lgl}
\alias{layout.lgl}
\concept{Graph layout}
\title{Large Graph Layout}
\description{
  A layout generator for larger graphs.
}
\usage{
layout.lgl (graph, maxiter=150, maxdelta=vcount(graph),
            area=vcount(graph)^2, coolexp=1.5,
            repulserad=area * vcount(graph),
            cellsize=sqrt(sqrt(area)), root=NULL)
}
\arguments{
  \item{graph}{The input graph}
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
\details{
  \code{layout.lgl} is for large connected graphs, it is similar to the
  layout generator of the Large Graph Layout software
  (\url{http://lgl.sourceforge.net/}).
}
\value{
  A numeric matrix with two columns and as many rows as vertices.
}
%\references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
%\seealso{}
% \examples{}
\keyword{graphs}
