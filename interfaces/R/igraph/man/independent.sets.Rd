\name{independent.vertex.sets}
\alias{independent.vertex.sets}
\alias{largest.independent.vertex.sets}
\alias{maximal.independent.vertex.sets}
\alias{independence.number}
\concept{Independent vertex set}
\title{Independent vertex sets}
\description{A vertex set is called independent if there no edges
  between any two vertices in it. These functions find independent
  vertex sets in undirected graphs}
\usage{
independent.vertex.sets(graph, min=NULL, max=NULL)
largest.independent.vertex.sets(graph)
maximal.independent.vertex.sets(graph)
independence.number(graph)
}
\arguments{
  \item{graph}{The input graph, directed graphs are considered as
    undirected, loop edges and multiple edges are ignored.}
  \item{min}{Numeric constant, limit for the minimum size of the
    independent vertex sets to find. \code{NULL} means no limit.}
  \item{max}{Numeric constant, limit for the maximum size of
    the independent vertex sets to find. \code{NULL} means no limit.
  }
}
\details{
  \code{independent.vertex.sets} finds all independent vertex sets in
  the network, obeying the size limitations given in the \code{min} and
  \code{max} arguments.

  \code{largest.independent.vertex.sets} finds the largest independent
  vertex sets in the graph. An independent vertex set is largest if
  there is no independent vertex set with more vertices.

  \code{maximal.independent.vertex.sets} finds the maximal independent
  vertex sets in the graph. An independent vertex set is maximal if it
  cannot be extended to a larger independent vertex set. The largest
  independent vertex sets are maximal, but the opposite is not always
  true.

  \code{independece.number} calculate the size of the largest
  independent vertex set(s).

  These functions use the algorithm described by Tsukiyama et al., see
  reference below.
}
\value{
  \code{independent.vertex.sets}, \code{largest.independent.vertex.sets}
  and \code{maximal.independent.vertex.sets} return a list containing
  numeric vertex ids, each list element is an independent vertex set.

  \code{independence.number} returns an integer constant.
}
\references{
  S. Tsukiyama, M. Ide, H. Ariyoshi and I. Shirawaka. A new algorithm
  for generating all the maximal independent sets. \emph{SIAM J Computing},
  6:505--517, 1977.
}
\author{Tamas Nepusz \email{ntamas@rmki.kfki.hu} ported it from the Very
  Nauty Graph Library by Keith Briggs \email{google@for.it} and Gabor Csardi
  \email{csardi@rmki.kfki.hu} wrote the R interface and this manual page. }
\seealso{\code{\link{cliques}}}
\examples{
# A quite dense graph
g <- erdos.renyi.game(100, 0.8)
independence.number(g)
independent.vertex.sets(g, min=independence.number(g))
largest.independent.vertex.sets(g)
# Empty graph
subgraph(g, largest.independent.vertex.sets(g)[[1]])

length(maximal.independent.vertex.sets(g))
}
\keyword{graphs}
