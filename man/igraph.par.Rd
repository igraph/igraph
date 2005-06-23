\name{igraph.par}
\alias{igraph.par}
\title{Igraph parameters}
\description{This functions can be used to query and set the parameters
  of the igraph package.}
\usage{
igraph.par(parid, parvalue = NULL)
}
\arguments{
  \item{parid}{String. The name of the parameter. See details below for
    a list of parameters.}
  \item{parvalue}{The new value of the parameter.}
}
\details{The package parameters:
  \itemize{default.type} The default storage type. This is
  \code{adjacencylist} right now.
}
\value{
  The value of the parameter for a query (\code{parvalue} argument not
  present. 

  TRUE when a parameter was set (\code{parvalue} argument present).
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso{}
\examples{
g <- graph.ring(10)
igraph.type(g)
igraph.par("default.type")
}
\keyword{graphs}
