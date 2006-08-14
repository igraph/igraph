\name{measure.dynamics}
\alias{measure.dynamics}
\alias{measure.dynamics.idage}
\alias{measure.dynamics.id}
\alias{measure.dynamics.d.d}
\alias{measure.dynamics.citedcat.id.age}
\alias{measure.dynamics.citingcat.id.age}
\title{Measuring the driving force in evolving networks}
\description{These functions assume a simple evolving network model and
  measure the functional form of a so-called \emph{attractiveness
    function} governing the evolution of the network.
}
\usage{
measure.dynamics.idage (graph, start.vertex = 0, agebins = 300, iterations = 5, 
    significance = 0, estind = NULL, estage = NULL, number = FALSE, 
    time.window = NULL)
measure.dynamics.id(graph, start.vertex = 0, iterations = 5, significance = 0, 
    estind = NULL, estage = NULL, number = FALSE, time.window = NULL) 
measure.dynamics.d.d(graph, vtime, etime, iterations = 5,
    sd = TRUE, no = FALSE) 
measure.dynamics.citedcat.id.age(graph, categories, agebins = 300,
    iterations = 5, significance = 0, 
    number = FALSE, norm = c(1, 1, 1)) 
measure.dynamics.citingcat.id.age(graph, categories, agebins = 300,
    iterations = 5, significance = 0,  
    number = FALSE, norm = c(1, 1, 1))
}
\arguments{
  \item{graph}{The graph of which the evolution is quantified. It is
    assumed that the vertices were added in increasing order of vertex
    id.}
  \item{start.vertex}{}
  \item{agebins}{Numeric constant, the number of bins to use for
    measuring aging.}
  \item{iterations}{Numeric constant, number of iterations to perform
    while calculating the attractiveness and the total attractiveness
    function.}
  \item{significance}{}
  \item{estind}{}
  \item{estage}{}
  \item{number}{}
  \item{time.window}{}
  \item{vtime}{}
  \item{etime}{}
  \item{sd}{}
  \item{no}{}
  \item{categories}{}
  \item{norm}{}  
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

  