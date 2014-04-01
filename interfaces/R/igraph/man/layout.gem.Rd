\name{layout.gem}
\alias{layout.gem}
\concept{Graph layout}
\title{The GEM layout algorithm}
\description{Place vertices on the plane using the GEM force-directed
  layout algorithm.}
\usage{
layout.gem (graph, coords=NULL, maxiter=40*vcount(graph)^2,
            temp.max=vcount(graph), temp.min=1/10,
            temp.init=sqrt(vcount(graph)))
}
\arguments{
  \item{graph}{The input graph. Edge directions are ignored.}
  \item{coords}{If not \code{NULL}, then the starting coordinates
    should be given here, in a two or three column matrix, depending on
    the \code{dim} argument.}
  \item{maxiter}{The maximum number of iterations to
    perform. Updating a single vertex counts as an iteration.
    A reasonable default is 40 * n * n, where n is the number of 
    vertices. The original paper suggests 4 * n * n, but this 
    usually only works if the other parameters are set up carefully.}
  \item{temp.max}{The maximum allowed local temperature. A reasonable
    default is the number of vertices.}
  \item{temp.min}{The global temperature at which the algorithm
    terminates (even before reaching \code{maxiter} iterations). A
    reasonable default is 1/10.}
  \item{temp.init}{Initial local temperature of all vertices. A
    reasonable default is the square root of the number of
    vertices.}
}
\details{
  See the referenced paper below for the details of the algorithm.
}
\value{
  A numeric matrix with two columns, and as many rows as the number of
  vertices.
}
\references{
  Arne Frick, Andreas Ludwig,  Heiko Mehldau: A Fast Adaptive Layout
  Algorithm for Undirected Graphs, \emph{Proc. Graph Drawing 1994}, LNCS
  894, pp. 388-403, 1995.
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{layout.fruchterman.reingold}},
  \code{\link{plot.igraph}}, \code{\link{tkplot}}}
\examples{
set.seed(42)
g <- graph.ring(10)
plot(g, layout=layout.gem)
}
\keyword{graphs}
