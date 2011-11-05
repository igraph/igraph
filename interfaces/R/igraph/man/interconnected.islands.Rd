\name{interconnected.islands}
\alias{interconnected.islands.game}
\concept{Random graph model}
\title{A graph with subgraphs that are each a random graph.}
\description{
  Create a number of Erdos-Renyi random graphs with identical
  parameters, and connect them with the specified number of edges.
}
\usage{
interconnected.islands.game (islands.n, islands.size, islands.pin,
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
\author{Samuel Thiriot \email{samuel.thiriot at res-ear dot ch}}
\seealso{\code{\link{erdos.renyi.game}}}
\examples{
g <- interconnected.islands.game(3, 10, 5/10, 1)
oc <- optimal.community(g)
oc
}
\keyword{graphs}

