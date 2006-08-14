\name{igraph-parameters}
\alias{igraph.par}
\title{Parameters for the igraph package}
\description{igraph has some parameters which (usually) affect the
  behavior of many functions. These can be set for the whole session via
  \code{igraph.par}.
}
\usage{
igraph.par(parid, parvalue = NULL)
}
\arguments{
  \item{parid}{The name of the parameter. See the currently used
    parameters below.}
  \item{parvalue}{The new value of the parameter. If \code{NULL} then
    the current value of the parameter is listed.}
}
\details{
  The currently used parameters in alphabetical order:
  \itemize{
    \item{print.edge.attributes}{Logical constant, whether to print edge
      attributes when printing graphs. Defaults to \code{FALSE}.}
    \item{print.graph.attributes}{Logical constant, whether to print
      graph attributes when printing graphs. Defaults to \code{FALSE}.}
    \item{print.vertex.attriubutes}{Logical constant, whether to print
      vertex attributes when printing graphs. Defaults to \code{FALSE}.}
    \item{verbose}{Logical constant, whether igraph functions should
      talk more than minimal. Eg. if \code{TRUE} thne some functions
      will use progress bars while computing. Defaults to \code{FALSE}.}
  }
}
\value{
  If \code{parvalue} is \code{NULL} then the current value of the
  parameter is returned. Otherwise the new value of the parameter is
  returned invisibly.
}
% \references{}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso{}
\examples{
igraph.par("verbose", FALSE)
layout.kamada.kawai( graph.ring(10) )
igraph.par("verbose", TRUE)
}
\keyword{graphs}
