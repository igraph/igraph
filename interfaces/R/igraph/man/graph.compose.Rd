\name{graph.compose}
\alias{graph.compose}
\alias{\%c\%}
\concept{Graph operators}
\title{Compose two graphs as binary relations}
\description{Relational composition of two graph.}
\usage{
graph.compose(g1, g2)
}
\arguments{
  \item{g1}{The first input graph.}
  \item{g2}{The second input graph.}
}
\details{
  \code{graph.compose} creates the relational composition of two
  graphs. The new graph will contain an (a,b) edge only if there is a
  vertex c, such that edge (a,c) is included in the first graph and
  (c,b) is included in the second graph. The corresponding operator is
  \%c\%.

  The function gives an error if one of the input graphs is directed and
  the other is undirected.

  Note that the function may generate multigraphs, if there are more
  than one way to find edges (a,b) in g1 and (b,c) in g2 for an edge
  (a,c) in the result. See \code{\link{simplify}} if you want to get rid
  of the multiple edges.

  The function may create loop edges, if edges (a,b) and (b,a) are
  present in g1 and g2, respectively, then (a,a) is included in the
  result. See \code{\link{simplify}} if you want to get rid of the
  self-loops.
}
\value{
  A new graph object.
}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
% \seealso{}
\examples{
g1 <- graph.ring(10)
g2 <- graph.star(10, mode="undirected")
gc <- graph.compose(g1, g2)
str(gc)
str(simplify(gc))
}
\keyword{graphs}
