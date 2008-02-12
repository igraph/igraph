\name{subgraph}
\alias{subgraph}
\concept{Subgraph}
\title{Subgraph of a graph}
\description{\code{subgraph} creates a subgraph of a graph, containing
  only the specified vertices and all the edges among them.}
\usage{
subgraph(graph, v)
}
\arguments{
  \item{graph}{The original graph.}
  \item{v}{Numeric vector, the vertices of the original graph which
    will form the subgraph.}
}
\details{
  The ids of the vertices will change in the subgraph of course since
  these are always consecutive.
}
\value{A new graph object.}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso
\examples{
g <- graph.ring(10)
g2 <- subgraph(g, 1:7)
}
\keyword{graphs}
