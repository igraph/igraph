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
as.undirected(graph, mode = c("collapse", "each"),
       edge.attr.comb = igraph.par("edge.attr.comb"))
}
\arguments{
  \item{graph}{The graph to convert.}
  \item{mode}{Character constant, defines the conversion algorithm. For
    \code{as.directed} it can be \code{mutual} or \code{arbitrary}. For
    \code{as.undirected} it can be \code{each} or \code{collapse}. See
    details below.}
  \item{edge.attr.comb}{Specifies what to do with edge attributes, if
    \code{mode="collapse"}. In this case many edges might be mapped to a
    single one in the new graph, and their attributes are
    combined. Please see \code{\link{attribute.combination}} for details
    on this.}
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

# Combining edge attributes
g3 <- graph.ring(10, directed=TRUE, mutual=TRUE)
E(g3)$weight <- seq_len(ecount(g3))
ug3 <- as.undirected(g3)
print(ug3, e=TRUE)
\dontrun{
  x11(width=10, height=5)
  layout(rbind(1:2))
  plot( g3, layout=layout.circle, edge.label=E(g3)$weight)
  plot(ug3, layout=layout.circle, edge.label=E(ug3)$weight)
}
}
\keyword{graphs}
