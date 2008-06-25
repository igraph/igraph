\name{is.mutual}
\alias{is.mutual}
\concept{Mutual edges}
\concept{Reciprocity}
\title{Find mutual edges in a directed graph}
\description{This function checks the reciproc pair of the supplied edges.}
\usage{
is.mutual(graph, es = E(graph))
}
\arguments{
  \item{graph}{The input graph.}
  \item{es}{Edge sequence, the edges that will be probed. By default is
    includes all edges in the order of their ids.}
}
\details{
  In a directed graph an (A,B) edge is mutual if the graph also includes
  a (B,A) directed edge.

  Note that multi-graphs are not handled properly, i.e. if the graph
  contains two copies of (A,B) and one copy of (B,A), then these three
  edges are considered to be mutual.

  Undirected graphs contain only mutual edges by definition.
}
\value{
  A logical vector of the same length as the number of edges supplied.
}
%\references{}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{ \code{\link{reciprocity}}, \code{\link{dyad.census}} if you
  just want some statistics about mutual edges.}
\examples{
g <- erdos.renyi.game(10,50,type="gnm",directed=TRUE)
reciprocity(g)
dyad.census(g)
is.mutual(g)
sum(is.mutual(g))/2 == dyad.census(g)$mut
}
\keyword{graphs}
