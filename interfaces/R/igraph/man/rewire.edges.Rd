\name{rewire.edges}
\alias{rewire.edges}
\title{Rewires the endpoints of the edges of a graph randomly}
\description{This function rewires the endpoints of the edges with a
  constant probability uniformly randomly to a new vertex in a graph. } 
\usage{
rewire.edges(graph, prob, loops=FALSE)
}
\arguments{
  \item{graph}{The input graph}
  \item{prob}{The rewiring probability, a real number between zero and
    one.}
  \item{loops}{Logical scalar, whether loop edges are allowed in the
    rewired graph.}
}
\details{
  Note that this function might create graphs with multiple and/or loop
  edges.
}
\value{A new graph object.}
%\references{}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu}}
%\seealso{}
\examples{
# Some random shortcuts shorten the distances on a lattice
g <- graph.lattice( length=100, dim=1, nei=5 )
average.path.length(g)
g <- rewire.edges( g, prob=0.05 )
average.path.length(g)
}
\keyword{graphs}
