\name{measure.dynamics}
\alias{measure.dynamics}
\alias{measure.dynamics.idage}
\alias{measure.dynamics.id}
\alias{measure.dynamics.d.d}
\alias{measure.dynamics.citedcat.id.age}
\alias{measure.dynamics.citingcat.id.age}
\alias{measure.dynamics.lastcit}
\alias{measure.dynamics.age}
\alias{measure.dynamics.citedcat}
\alias{measure.dynamics.citingcat.citedcat}
\title{Measuring the driving force in evolving networks}
\description{These functions assume a simple evolving network model and
  measure the functional form of a so-called \emph{attractiveness
    function} governing the evolution of the network.
}
\usage{
measure.dynamics.idage (graph, agebins = 300, iterations = 5, 
    error=TRUE, time.window = NULL, number = FALSE, cites=FALSE,
    norm.method="old")
measure.dynamics.id(graph, iterations = 5, error=TRUE, 
    time.window = NULL, number = FALSE, cites=FALSE, norm.method="old",
    debug=FALSE, debugdeg=0, which=2) 
measure.dynamics.d.d(graph, vtime, etime, iterations = 5)
measure.dynamics.citedcat.id.age(graph, categories, agebins = 300,
    iterations = 5, norm = c(1, 1, 1)) 
measure.dynamics.citingcat.id.age(graph, categories, agebins = 300,
    iterations = 5, norm = c(1, 1, 1))
measure.dynamics.lastcit(graph, agebins, iterations=5,
    norm.method="old", number=FALSE)
measure.dynamics.age(graph, agebins, iterations=5, norm.method="old",
    number=FALSE)
measure.dynamics.citedcat(graph, categories, iterations=5,
    number=FALSE, norm.method="old")
measure.dynamics.citingcat.citedcat(graph, categories, iterations=5,
    number=FALSE, norm.method="old", norm=c(1,1))

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
  \item{time.window}{}
  \item{vtime}{}
  \item{etime}{}
  \item{categories}{}
  \item{norm}{}
  \item{number}{}
  \item{norm.method}{}
  \item{error}{}
  \item{cites}{}
}
\details{
  The functions should be considered as experimental, so no
  detailed documentation yet. Sorry.
}
\value{
  TODO
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
% seealso{}
% \examples{}
\keyword{graphs}

  