\name{erdos.renyi.game}
\alias{erdos.renyi.game}
\alias{random.graph.game}
\title{Generate random graphs according to the Erdos-Renyi model}
\description{This model is very simple, every possible edge is created
  with the same constant probability. This is only one version of the ER
  model, the other -- adding the specified number of edges randomly to
  ampty graph -- is not implemented yet.}
\usage{
erdos.renyi.game(n, p, directed = FALSE, loops = FALSE, \dots)
}
\arguments{
  \item{n}{The number of vertices in the graph.}
  \item{p}{The probability for drawing an edge between two arbitrary
    vertices.}
  \item{directed}{Logical, whether the graph will be directed, defaults
    to FALSE.}
  \item{loops}{Logical, whether to add loop edges, defaults to FALSE.}
}
\details{
  The additional arguments are added as graph arguments.

  \code{random.graph.game} is an alias to this function.
}
\value{A graph object.}
\references{ Erdos, P. and Renyi, A., On random graphs, Publicationes
  Mathematicae 6, 290­297 (1959).
}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{barabasi.game}}}
\examples{
g <- erdos.renyi.game(1000, 1/1000)
degree.distribution(g)
}
\keyword{graphs}
