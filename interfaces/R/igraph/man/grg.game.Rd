\name{grg.game}
\alias{grg.game}
\concept{Geometric random graph}
\title{Geometric random graphs}
\description{Generate a random graph based on the distance of random
  point on a unit square}
\usage{
grg.game(nodes, radius, torus=FALSE, coords=FALSE)
}
\arguments{
  \item{nodes}{The number of vertices in the graph.}
  \item{radius}{The radius within which the vertices will be connected
    by an edge.}
  \item{torus}{Logical constant, whether to use a torus instead of a
    square.}
  \item{coords}{Logical scalar, whether to add the positions of the
    vertices as vertex attributes called \sQuote{\code{x}} and
    \sQuote{\code{y}}.}
}
\details{
  First a number of points are dropped on a unit square, these points
  correspond to the vertices of the graph to create. Two points will be
  connected with an undirected edge if they are closer to each other in
  Euclidean norm than a given radius. If the \code{torus} argument is
  \code{TRUE} then a unit area torus is used instead of a square.
}
\value{
  A graph object. If \code{coords} is \code{TRUE} then with vertex
  attributes \sQuote{\code{x}} and \sQuote{\code{y}}.
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}, first version was
  written by Keith Briggs \email{keith.briggs@bt.com}}
\seealso{\code{\link{random.graph.game}}}
\examples{
g <- grg.game(1000, 0.05, torus=FALSE)
g2 <- grg.game(1000, 0.05, torus=TRUE)
}
\keyword{graphs}
