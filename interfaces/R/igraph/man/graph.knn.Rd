\name{graph.knn}
\alias{graph.knn}
\title{Average nearest neighbor degree}
\description{Calculate the average nearest neighbor degree of the given
  vertices and the same quantity in the function of vertex degree}
\usage{
graph.knn(graph, vids=V(graph), weights=NULL)
}
\arguments{
  \item{graph}{The input graph. It can be directed, but it will be
    treated as undirected, i.e. the direction of the edges is ignored.}
  \item{vids}{The vertices for which the calculation is
    performed. Normally it includes all vertices. Note, that if not all
    vertices are given here, then both \sQuote{\code{knn}} and
    \sQuote{\code{knnk}} will be calculated based on the given vertices
    only.}
  \item{weights}{Weight vector. If the graph has a \code{weight} edge
    attribute, then this is used by default. If this argument is given,
    then vertex strength (see \code{\link{graph.strength}}) is used
    instead of vertex degree. But note that \code{knnk} is still given
    in the function of the normal vertex degree. }
}
\details{
  Note that for zero degree vertices the answer in \sQuote{\code{knn}}
  is \code{NaN} (zero divided by zero), the same is true for
  \sQuote{\code{knnk}} if a given degree never appears in the network.
}
\value{
  A list with two members:
  \item{knn}{A numeric vector giving the average nearest neighbor
    degree for all vertices in \code{vids}.}
  \item{knnk}{A numeric vector, its length is the maximum (total)
    vertex degree in the graph. The first element is the
    average nearest neighbor degree of vertices with degree one, etc.
  }
}
\references{
  Alain Barrat, Marc Barthelemy, Romualdo Pastor-Satorras, Alessandro
  Vespignani: The architecture of complex weighted networks,
  Proc. Natl. Acad. Sci. USA 101, 3747 (2004)
}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\examples{
# Some trivial ones
g <- graph.ring(10)
graph.knn(g)
g2 <- graph.star(10)
graph.knn(g2)

# A scale-free one, try to plot 'knnk'
g3 <- simplify(ba.game(1000, m=5))
graph.knn(g3)

# A random graph
g4 <- random.graph.game(1000, p=5/1000)
graph.knn(g4)

# A weighted graph
g5 <- graph.star(10)
E(g5)$weight <- seq(ecount(g5))
graph.knn(g5)
}
\keyword{graphs}
