\name{line.graph}
\alias{line.graph}
\concept{Line graph}
\title{Line graph of a graph}
\description{This function calculates the line graph of another graph.}
\usage{
line.graph(graph)
}
\arguments{
  \item{graph}{The input graph, it can be directed or undirected.}
}
\details{
 The line graph \code{L(G)} of a \code{G} undirected graph is defined as
 follows. \code{L(G)} has one vertex for each edge in \code{G} and two
 vertices in \code{L(G)} are connected  by an edge if their
 corresponding edges share an end point.

 The line graph \code{L(G)} of a \code{G} directed graph is slightly
 different, \code{L(G)} has one vertex for each edge in \code{G} and two
 vertices in \code{L(G)} are connected by a directed edge if the target
 of the first vertex's corresponding edge is the same as the source of
 the second vertex's corresponding edge.
}
\value{
  A new graph object.
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}, the first version of the C
  code was written by Vincent Matossian.}
% \seealso{}
\examples{
# generate the first De-Bruijn graphs
g <- graph.full(2, directed=TRUE, loops=TRUE)
line.graph(g)
line.graph(line.graph(g))
line.graph(line.graph(line.graph(g)))
}
\keyword{graphs}
