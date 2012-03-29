\name{traits}
\alias{callaway.traits.game}
\alias{establishment.game}
\title{Graph generation based on different vertex types}
\description{These functions implement evolving network models based on
  different vertex types.
}
\usage{
callaway.traits.game (nodes, types, edge.per.step = 1, type.dist = rep(1, 
    types), pref.matrix = matrix(1, types, types), directed = FALSE) 
establishment.game(nodes, types, k = 1, type.dist = rep(1, types),
    pref.matrix = matrix(1, types, types), directed = FALSE) 
}
\arguments{
  \item{nodes}{The number of vertices in the graph.}
  \item{types}{The number of different vertex types.}
  \item{edge.per.step}{The number of edges to add to the graph per time step.}
  \item{type.dist}{The distribution of the vertex types. This is assumed
    to be stationary in time.}
  \item{pref.matrix}{A matrix giving the preferences of the given vertex
    types. These should be probabilities, ie. numbers between zero and
    one.}
  \item{directed}{Logical constant, whether to generate directed
    graphs.}
  \item{k}{The number of trials per time step, see details below.}
}
\details{
  For \code{callaway.traits.game} the simulation goes like this: in each
  discrete time step a new vertex is added to the graph. The type of
  this vertex is generated based on \code{type.dist}. Then two vertices are
  selected uniformly randomly from the graph. The probability that they
  will be connected depends on the types of these vertices and is taken
  from \code{pref.matrix}. Then another two vertices are selected and
  this is repeated \code{edges.per.step} times in each time step.

  For \code{establishment.game} the simulation goes like this: a single
  vertex is added at each time step. This new vertex tries to connect to
  \code{k} vertices in the graph. The probability that such a connection is
  realized depends on the types of the vertices involved and is taken
  from \code{pref.matrix}.
}
\value{
  A new graph object.
}
% \references{}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso{}
\examples{
# two types of vertices, they like only themselves
g1 <- callaway.traits.game(1000, 2, pref.matrix=matrix( c(1,0,0,1), ncol=2))
g2 <- establishment.game(1000, 2, k=2, pref.matrix=matrix( c(1,0,0,1), ncol=2))
}
\keyword{graphs}
