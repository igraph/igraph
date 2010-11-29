\name{is.weighted}
\alias{is.weighted}
\concept{Weighted graphs}
\title{Weighted graphs}
\description{In weighted graphs, a real number is assigned to each
  (directed or undirected) edge.}
\usage{
is.weighted(graph)
}
\arguments{
  \item{graph}{The input graph.}
}
\details{
  In igraph edge weights are represented via an edge attribute, called
  \sQuote{weight}. The \code{is.weighted} function only checks that such
  an attribute exists. (It does not even checks that it is a numeric edge
  attribute.)

  Edge weights are used for different purposes by the different
  functions. E.g. shortest path functions use it as the cost of the
  path; community finding methods use it as the strength of the
  relationship between two vertices, etc. Check the manual pages of the
  functions working with weighted graphs for details.
}
\value{A logical scalar.}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\examples{
g <- graph.ring(10)
get.shortest.paths(g, 8, 2)
E(g)$weight <- seq_len(ecount(g))
get.shortest.paths(g, 8, 2)
}
\keyword{graphs}
