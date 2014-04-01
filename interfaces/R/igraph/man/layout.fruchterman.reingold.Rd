\name{layout.fruchterman.reingold}
\alias{layout.fruchterman.reingold}
\concept{Graph layout}
\title{The Fruchterman-Reingold layout algorithm}
\description{Place vertices on the plane using the force-directed
  layout algorithm by Fruchterman and Reingold.}
\usage{
layout.fruchterman.reingold (graph, coords = NULL, dim = 2,
    niter = 500, start.temp = sqrt(vcount(graph)), grid=c("auto",
      "grid", "nogrid"), weights = NULL, minx = NULL, maxx = NULL,
    miny = NULL, maxy = NULL, minz = NULL, maxz = NULL,
    coolexp, maxdelta, area, repulserad)
}
\arguments{
  \item{graph}{The graph to lay out. Edge directions are ignored.}
  \item{coords}{Optional starting positions for the vertices. If this
    argument is not \code{NULL} then it should be an appropriate
    matrix of starting coordinates.}
  \item{dim}{Integer scalar, 2 or 3, the dimension of the layout.
    Two dimensional layouts are places on a plane, three dimensional
    ones in the 3d space.}
  \item{niter}{Integer scalar, the number of iterations to perform.}
  \item{start.temp}{Real scalar, the start temperature. This is the
    maximum amount of movement alloved along one axis, within one step,
    for a vertex. Currently it is decreased linearly to zero during
    the iteration.}
  \item{grid}{Character scalar, whether to use the faster, but less
    accurate grid based implementation of the algorithm. By default
    (\dQuote{auto}), the grid-based implementation is used if the graph
    has more than one thousand vertices.}
  \item{weights}{A vector giving edge weights. The \code{weight} edge
    attribute is used by default, if present. If weights are given, then
    the attraction along the edges will be multiplied
    by the given edge weights.}
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
  \item{coolexp,maxdelta,area,repulserad}{These arguments are not
    supported from igraph version 0.8.0 and are ignored (with a warning).}
}
\details{
  See the referenced paper below for the details of the algorithm.

  This function was rewritten from scratch in igraph version 0.8.0.
}
\value{
  A two- or three-column matrix, each row giving the coordinates of a
  vertex, according to the ids of the vertex ids.
}
\references{
  Fruchterman, T.M.J. and Reingold, E.M. (1991). Graph Drawing by
  Force-directed Placement. \emph{Software - Practice and Experience},
  21(11):1129-1164. 
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{layout.drl}}, \code{\link{layout.kamada.kawai}} for
  other layout algorithms.}
\examples{
# Fixing ego
g <- ba.game(20, m=2)
minC <- rep(-Inf, vcount(g))
maxC <- rep(Inf, vcount(g))
minC[1] <- maxC[1] <- 0
co <- layout.fruchterman.reingold(g, minx=minC, maxx=maxC,
                                  miny=minC, maxy=maxC)
co[1,]
plot(g, layout=co, vertex.size=30, edge.arrow.size=0.2,
     vertex.label=c("ego", rep("", vcount(g)-1)), rescale=FALSE,
     xlim=range(co[,1]), ylim=range(co[,2]), vertex.label.dist=0,
     vertex.label.color="red")
axis(1)
axis(2)
}
\keyword{graphs}
