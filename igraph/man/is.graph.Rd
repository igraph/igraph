\name{is.igraph}
\alias{is.igraph}
\title{Is this object a graph?}
\description{\code{is.graph} makes its decision based on the class
  attribute of the object.}
\usage{
is.igraph(graph)
}
\arguments{
  \item{graph}{An R object.}
}
% \details{}
\value{
  A logical constant, \code{TRUE} if argument \code{graph} is a graph
  object.
}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso{}
\examples{
g <- graph.ring(10)
is.igraph(g)
is.igraph(numeric(10))
}
\keyword{graphs}
