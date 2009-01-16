\name{girth}
\alias{girth}
\concept{Girth}
\title{Girth of a graph}
\description{The girth of a graph is the length of the shortest circle
  in it.}
\usage{
girth(graph, circle=TRUE)
}
\arguments{
  \item{graph}{The input graph. It may be directed, but the algorithm
    searches for undirected circles anyway.}
  \item{circle}{Logical scalar, whether to return the shortest circle
    itself.}
}
\details{
  The current implementation works for undirected graphs only, 
  directed graphs are treated as undirected graphs. Loop edges and
  multiple edges are ignored.
  If the graph is a forest (ie. acyclic), then zero is returned.

  This implementation is based on Alon Itai and Michael Rodeh:
  Finding a minimum circuit in a graph
  \emph{Proceedings of the ninth annual ACM symposium on Theory of
    computing}, 1-10, 1977. The first implementation of this
  function was done by Keith Briggs, thanks Keith.
}
\value{
  A named list with two components:
  \item{girth}{Integer constant, the girth of the graph, or 0 if the
    graph is acyclic.}
  \item{circle}{Numeric vector with the vertex ids in the shortest
    circle.}
}
\references{
  Alon Itai and Michael Rodeh: Finding a minimum circuit in a graph 
  \emph{Proceedings of the ninth annual ACM symposium on Theory of 
    computing}, 1-10, 1977
}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\examples{
# No circle in a tree
g <- graph.tree(1000, 3)
girth(g)

# The worst case running time is for a ring
g <- graph.ring(100)
girth(g)

# What about a random graph?
g <- erdos.renyi.game(1000, 1/1000)
girth(g)
}
\keyword{graphs}

