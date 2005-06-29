\name{cocitation}
\alias{cocitation}
\title{Cocitation coupling}
\description{Two vertices are cocited if there is another vertex citing
  both of them. This function siply counts how many types two vertices
  are cocited.}
\usage{
cocitation(graph, v = 1:vcount(graph))
}
\arguments{
  \item{graph}{The graph object to analyze}
  \item{v}{Numeric vector, the vertex ids for which the cocitation
    counts we want to calculate. The default is all vertices.}
}
\details{
  The function calculates the cocitation counts for the vertices in the
  \code{v} argument and all vertices in the graph.

  Calculating the cocitation for only one vertex costs the same amount
  of computation as for all vertices. This might change in the future.
}
\value{
  A numeric matrix with \code{length(v)} lines and \code{vcount(graph)}
  columns. Element \code{(i,j)} contains the number of the times a
  vertex cited both \code{i} and \code{j}.
}
% \references{}
\author{Gabor Csardi \email{csrdi@rmki.kfki.hu}}
% \seealso{}
\examples{
g <- graph.ring(10)
cocitation(g)
}
\keyword{graphs}

