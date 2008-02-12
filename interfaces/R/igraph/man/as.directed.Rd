\name{as.directed}
\alias{as.directed}
\alias{as.undirected}
\concept{Directed graph}
\concept{Undirected graph}
\title{Convert between directed and undirected graphs}
\description{\code{as.directed} converts an undirected graph to
  directed, \code{as.undirected} is the opposite, it converts a directed
  graph to undirected.}
\usage{
as.directed(graph, mode = c("mutual", "arbitrary")) 
as.undirected(graph, mode = c("collapse", "each")) 
}
\arguments{
  \item{graph}{The graph to convert.}
  \item{mode}{Character constant, defines the conversion algorithm. For
    \code{as.directed} it can be \code{mutual} or \code{arbitrary}. For
    \code{as.undirected} it can be \code{each} or \code{collapse}. See
    details below.}
}
\details{
  Conversion algorithms for \code{as.directed}:
  \itemize{
    \item{arbitrary}{The number of
      edges in the graph stays the same, an arbitrarily directed edge is
      created for each undirected edge.}
    \item{mutual}{Two directed edges are created for each undirected
      edge, one in each direction.}
  }

  Conversion algorithms for \code{as.undirected}:
  \itemize{
    \item{each}{The number of edges remains constant, an undirected edge
      is created for each directed one, this version might create graphs
      with multiple edges.}
    \item{collapse}{One undirected edge will be created for each pair of
      vertices which are connected with at least one directed edge, no
      multiple edges will be created.}
  }
}
\value{
  A new graph object.
}
%\references{}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{simplify}} for removing multiple and/or loop edges
  from a graph.}
\examples{
g <- graph.ring(10)
as.directed(g, "mutual")
g2 <- graph.star(10)
as.undirected(g)
}
\keyword{graphs}
