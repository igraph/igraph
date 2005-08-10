\name{aging.prefatt.game}
\alias{aging.prefatt.game}
\title{Generate an evolving random graph with preferential attachment
  and aging}
\description{This function creates a random graph by simulation its
  evolution. Each time a new vertex is added it creates a number of links
  to old vertices and the probability that an old vertex is cited depends on
  its in-degree (preferential attachment) and age.}
\usage{
aging.prefatt.game(n, m=1, aging.type="exponential", params=list(), \dots)
}
\arguments{
  \item{n}{The number of vertices in the graph.}
  \item{m}{The number of edges each new vertex creates (except the very
    first vertex.}
  \item{aging.type}{Character string, the type of the function giving
    probability that an old vertex is cited depending on its age. See
    details below.}
  \item{params}{Named list, this gives the parameters of the
    \dQuote{aging function} selected by the \code{aging.type}
    argument. See details below.}
  \item{\dots}{Additional arguments, these are passed to the graph
    constructor.}
}
\details{
  This is discrete time step model, in each time step a new vertex is
  added to the network. The new vertex cites a number (parameter
  \code{m}) of other vertices. The probability that a vertex is cited is
  proportional to the product of the in-degree of the node plus one and
  the so-called aging function (\code{aging.type} arguments).

  The \code{exponential} aging function decreases exponetially with age
  with exponent \code{aging.exp}, its only parameter.

  The \code{powerlaw} aging function decreases as a power law with age,
  the exponent is given by the \code{aging.exp} parameter.
}
\value{A new graph.}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{ba.game}}, \code{\link{erdos.renyi.game}}}
\examples{
g <- aging.prefatt.game(100, aging.type="exponential",
     params=list(aging.exp=1/100))
g2 <- aging.prefatt.game(100, aging.type="powerlaw",
     params=list(aging.exp=1))
}
\keyword{graphs}
