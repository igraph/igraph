\name{watts.strogatz.game}
\alias{watts.strogatz.game}
\concept{Small-world model}
\concept{Watts-strogatz model}
\title{The Watts-Strogatz small-world model}
\description{Generate a graph according to the Watts-Strogatz network
  model.}
\usage{
watts.strogatz.game(dim, size, nei, p)
}
\arguments{
  \item{dim}{Integer constant, the dimension of the starting lattice.}
  \item{size}{Integer constant, the size of the lattice along each
    dimension.}
  \item{nei}{Integer constant, the neighborhood within which the
    vertices of the lattice will be connected.}
  \item{p}{Real constant between zero and one, the rewiring probability. }
}
\details{
  First a lattice is created with the given \code{dim}, \code{size} and
  \code{nei} arguments. Then the edges of the lattice are rewired
  uniformly randomly with probability \code{p}.

  Note that this function might create graphs with loops and/or multiple
  edges. You can use \code{\link{simplify}} to get rid of these.
}
\value{A graph object.}
\references{
  Duncan J Watts and Steven H Strogatz: Collective dynamics of
  \sQuote{small world} networks, Nature 393, 440-442, 1998.
}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{graph.lattice}}, \code{\link{rewire.edges}}}
\examples{
g <- watts.strogatz.game(1, 100, 5, 0.05)
average.path.length(g)
transitivity(g, type="average")
}
\keyword{graphs}
