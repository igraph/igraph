\name{graph.compose}
\alias{graph.compose}
\alias{\%c\%}
\concept{Graph operators}
\title{Compose two graphs as binary relations}
\description{Relational composition of two graph.}
\usage{
graph.compose(g1, g2, byname = "auto")
}
\arguments{
  \item{g1}{The first input graph.}
  \item{g2}{The second input graph.}
  \item{byname}{A logical scalar, or the character scalar
    \code{auto}. Whether to perform the operation based on symbolic
    vertex names. If it is \code{auto}, that means \code{TRUE} if both
    graphs are named and \code{FALSE} otherwise. A warning is generated
    if \code{auto} and one graph, but not both graphs are named.}
}
\details{
  \code{graph.compose} creates the relational composition of two
  graphs. The new graph will contain an (a,b) edge only if there is a
  vertex c, such that edge (a,c) is included in the first graph and
  (c,b) is included in the second graph. The corresponding operator is
  \%c\%.

  The function gives an error if one of the input graphs is directed and
  the other is undirected.

  If the \code{byname} argument is \code{TRUE} (or \code{auto} and the
  graphs are all named), then the operation is performed based on
  symbolic vertex names. Otherwise numeric vertex ids are used.

  \code{graph.compose} keeps the attributes of both graphs. All
  graph, vertex and edge attributes are copied to the result. If an
  attribute is present in multiple graphs and would result a name clash,
  then this attribute is renamed by adding suffixes: _1, _2, etc.

  The \code{name} vertex attribute is treated specially if the operation
  is performed based on symbolic vertex names. In this case \code{name}
  must be present in both graphs, and it is not renamed in the result
  graph.

  Note that an edge in the result graph corresponds to two edges in the
  input, one in the first graph, one in the second. This mapping is not
  injective and several edges in the result might correspond to the same
  edge in the first (and/or the second) graph. The edge attributes in
  the result graph are updated accordingly.

  Also note that the function may generate multigraphs, if there are more
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
