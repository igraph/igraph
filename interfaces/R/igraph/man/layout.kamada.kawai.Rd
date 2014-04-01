\name{layout.kamada.kawai}
\alias{layout.kamada.kawai}
\concept{Graph layout}
\title{The Kamada-Kawai layout algorithm}
\description{Place the vertices on the plane, or in the 3d space,
  based on a phyisical model of springs.}
\usage{
layout.kamada.kawai (graph, coords = NULL, dim = 2,
    maxiter = 50 * vcount(graph), epsilon = 0, kkconst = vcount(graph),
    weights = NULL, minx = NULL, maxx = NULL, miny = NULL, maxy = NULL,
    minz = NULL, maxz = NULL, niter, sigma, initemp, coolexp)
}
\arguments{
  \item{graph}{The input graph. Edge directions are ignored.}
  \item{coords}{If not \code{NULL}, then the starting coordinates
    should be given here, in a two or three column matrix, depending on
    the \code{dim} argument.}
  \item{dim}{Integer scalar, 2 or 3, the dimension of the layout.
    Two dimensional layouts are places on a plane, three dimensional
    ones in the 3d space.}
  \item{maxiter}{The maximum number of iterations to perform. The
    algorithm might terminate earlier, see the \code{epsilon} argument.}
  \item{epsilon}{Numeric scalar, the algorithm terminates, if the
    maximal delta is less than this. (See the reference below for what
    delta means.) If you set this to zero, then the function always
    performs \code{maxiter} iterations.}
  \item{kkconst}{Numeric scalar, the Kamada-Kawai vertex attraction
    constant. Typical (and default) value is the number of vertices.}
  \item{weights}{Edge weights, larger values will result longer edges.}
  \item{minx}{If not \code{NULL}, then it must be a numeric vector
    that gives lower boundaries for the \sQuote{x} coordinates of the
    vertices. The length of the vector must match the number of
    vertices in the graph.}
  \item{maxx}{Similar to \code{minx}, but gives the upper boundaries.}
  \item{miny}{Similar to \code{minx}, but gives the lower boundaries
    of the \sQuote{y} coordinates.}
  \item{maxy}{Similar to \code{minx}, but gives the upper boundaries
    of the \sQuote{y} coordinates.}
  \item{minz}{Similar to \code{minx}, but gives the lower boundaries
    of the \sQuote{z} coordinates.}
  \item{maxz}{Similar to \code{minx}, but gives the upper boundaries
    of the \sQuote{z} coordinates.}
  \item{niter,sigma,initemp,coolexp}{These arguments are not
    supported from igraph version 0.8.0 and are ignored (with a warning).}
}
\details{
  See the referenced paper below for the details of the algorithm.
  
  This function was rewritten from scratch in igraph version 0.8.0 and
  it follows truthfully the original publication by Kamada and Kawai
  now.
}
\value{
  A numeric matrix with two (dim=2) or three (dim=3) columns, and as
  many rows as the number of vertices, the x, y and potentially z
  coordinates of the vertices.
}
\references{
  Kamada, T. and Kawai, S.: An Algorithm for Drawing General Undirected
  Graphs. \emph{Information Processing Letters}, 31/1, 7--15, 1989.
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{layout.drl}}, \code{\link{plot.igraph}},
  \code{\link{tkplot}}}
\examples{
g <- graph.ring(10)
E(g)$weight <- rep(1:2, length.out=ecount(g))
plot(g, layout=layout.kamada.kawai, edge.label=E(g)$weight)
}
\keyword{graphs}

