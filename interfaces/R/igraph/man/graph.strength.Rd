\name{graph.strength}
\alias{graph.strength}
\title{Strength or weighted vertex degree}
\description{Summing up the edge weights of the adjacent edges for each
  vertex.}
\usage{
graph.strength (graph, vids = V(graph), mode = c("all", "out", "in", "total"),
      loops = TRUE, weights = NULL)
}
\arguments{
  \item{graph}{The input graph.}
  \item{vids}{The vertices for which the strength will be calculated.}
  \item{mode}{Character string, \dQuote{out} for out-degree, \dQuote{in} for
    in-degree or \dQuote{all} for the sum of the two. For undirected
    graphs this argument is ignored.}
  \item{loops}{Logical; whether the loop edges are also counted.}
  \item{weights}{Weight vector. If the graph has a \code{weight} edge
    attribute, then this is used by default. If the graph does not have
    a \code{weight} edge attribute and this argument is \code{NULL},
    then a warning is given and \code{\link{degree}} is called.}
}
% \details{}
\value{
  A numeric vector giving the strength of the vertices.
}
\references{
  Alain Barrat, Marc Barthelemy, Romualdo Pastor-Satorras, Alessandro
  Vespignani: The architecture of complex weighted networks,
  Proc. Natl. Acad. Sci. USA 101, 3747 (2004)
}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu} }
\seealso{\code{\link{degree}} for the unweighted version.}
\examples{
g <- graph.star(10)
E(g)$weight <- seq(ecount(g))
graph.strength(g)
graph.strength(g, mode="out")
graph.strength(g, mode="in")

# No weights, a warning is given
g <- graph.ring(10)
graph.strength(g)
}
\keyword{graphs}
