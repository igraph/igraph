\name{cocitation}
\alias{cocitation}
\alias{bibcoupling}
\title{Cocitation coupling}
\description{Two vertices are cocited if there is another vertex citing
  both of them. \code{cocitation} siply counts how many types two vertices
  are cocited. The bibliographic coupling of two vertices is the number
  of other vertices they both cite, \code{bibcoupling} calculates this.
}
\usage{
cocitation(graph, v = 0:(vcount(graph)-1))
bibcoupling(graph, v = 0:(vcount(graph)-1))
}
\arguments{
  \item{graph}{The graph object to analyze}
  \item{v}{Numeric vector, the vertex ids for which the cocitation or
    bibliographic coupling values we want to calculate. The default is
    all vertices.} 
}
\details{
  \code{cocitation} calculates the cocitation counts for the vertices in the
  \code{v} argument and all vertices in the graph.

  \code{bibcoupling} calculates the bibliographic coupling for vertices
  in \code{v} and all vertices in the graph.

  Calculating the cocitation or bibliographic coupling for only one
  vertex costs the same amount of computation as for all vertices. This
  might change in the future. 
}
\value{
  A numeric matrix with \code{length(v)} lines and \code{vcount(graph)}
  columns. Element \code{(i,j)} contains the cocitation or bibliographic
  coupling for vertices \code{i} and \code{j}.
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso{}
\examples{
g <- graph.ring(10)
cocitation(g)
bibcoupling(g)
}
\keyword{graphs}

