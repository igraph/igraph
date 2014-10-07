\name{interconnected.islands}
\alias{interconnected.islands.game}
\alias{sample_islands}
\concept{Random graph model}
\title{A graph with subgraphs that are each a random graph.}
\description{
  Create a number of Erdos-Renyi random graphs with identical
  parameters, and connect them with the specified number of edges.
}
\usage{
sample_islands (islands.n, islands.size, islands.pin,
                             n.inter)
}
\arguments{
  \item{islands.n}{The number of islands in the graph.}
  \item{islands.size}{The size of islands in the graph.}
  \item{islands.pin}{The probability to create each possible edge into
    each island.}
  \item{n.inter}{The number of edges to create between two islands.}
}
% \details{}
\value{
  An igraph graph.
}
% \references{}
\author{Samuel Thiriot (\url{http://samuelthiriot.res-ear.ch/})}
\seealso{\code{\link{sample_gnp}}}
\examples{
g <- sample_islands(3, 10, 5/10, 1)
oc <- cluster_optimal(g)
oc
}
\keyword{graphs}

