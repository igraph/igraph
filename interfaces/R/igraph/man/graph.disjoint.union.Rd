\name{graph.disjoint.union}
\alias{graph.disjoint.union}
\alias{\%du\%}
\concept{Graph operators}
\title{Disjoint union of graphs}
\description{The union of two or more graphs are created. The graphs are
  assumed to have disjoint vertex sets.}
\usage{
graph.disjoint.union(\dots)
}
\arguments{
  \item{\dots}{Graph objects or lists of graph objects.}
}
\details{
  \code{graph.disjoint.union} creates a union of two or more disjoint
  graphs. Thus first the vertices in the second, third, etc. graphs are
  relabeled to have completely disjoint graphs. Then a simple union is
  created. This function can also be used via the \%du\% operator.  
  
  \code{graph.disjont.union} handles graph, vertex and edge attributes.
  In particular, it merges vertex and edge
  attributes using the basic \code{c()} function. For graphs that lack
  some vertex/edge attribute, the corresponding values in the new graph
  are set to \code{NA}. Graph attributes are simply copied to the
  result. If this would result a name clash, then they are renamed by
  adding suffixes: _1, _2, etc.

  Note that if both graphs have vertex names (ie. a \code{name} vertex
  attribute), then the concatenated vertex names might be non-unique
  in the result. A warning is given if this happens.
  
  An error is generated if some input graphs are directed and others are
  undirected.
}
\value{
  A new graph object.
}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
% \seealso{}
\examples{
## A star and a ring
g1 <- graph.star(10, mode="undirected")
V(g1)$name <- letters[1:10]
g2 <- graph.ring(10)
V(g2)$name <- letters[11:20]
str(g1 \%du\% g2)
}
\keyword{graphs}
