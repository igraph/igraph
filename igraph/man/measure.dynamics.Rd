\name{measure.dynamics}
\alias{measure.dynamics}
\alias{measure.dynamics.idage}
\title{Measuring the driving force in evolving networks}
\description{These functions assume a simple evolving network model and
  measure the functional form of a so-called \emph{attractiveness
    function} governing the evolution of the network.
}
\usage{
measure.dynamics.idage(graph, agebins=300, iterations=5, sd=FALSE)
}
\arguments{
  \item{graph}{The graph of which the evolution is quantified. It is
    assumed that the vertices were added in increasing order of vertex
    id.} 
  \item{agebins}{Numeric constant, the number of bins to use for
    measuring aging.}
  \item{iterations}{Numeric constant, number of iterations to perform
    while calculating the attractiveness and the total attractiveness
    function.}
  \item{sd}{Logical, should an estimation of the error also calculated.}
}
\details{
  The functions should be considered as experimental, so no
  detailed documentation yet. Sorry.
}
\value{
  A named list with components \code{akl} (the attractiveness function),
  \code{st} (total attractiveness versus time), \code{sd} (only if the
  \code{sd} arguments was \code{TRUE}, the estimation of the error for
  the \code{akl} values).
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
% seealso{}
\examples{
g <- ba.game(10000)
mes <- measure.dynamics.idage(g, 20, sd=TRUE)
}
\keyword{graphs}

  