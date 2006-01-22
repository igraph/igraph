\name{transitivity}
\alias{transitivity}
\title{Transitivity of a graph}
\description{Transitivity measures the probability that the adjacent
  vertices of a vertex are connected. This is sometimes also called the
  clustering coefficient. }
\usage{
transitivity(graph, type="undirected")
}
\arguments{
  \item{graph}{The graph to analyze.}
  \item{type}{The type of the transitivity to calculate. Possible
    values: \sQuote{\code{undirected}} calculates the global transitivity of an
    undirected graph (directed graphs are considered as undirected ones
    as well). This is simply the ratio of the triangles and the
    connected triples in the graph.
  }
}
%\details{}
\value{
  For \sQuote{\code{undirected}} a single number, or \code{NaN} if there are no
  connected triples in the graph.
}
\references{ Wasserman, S., and Faust, K. (1994). \emph{Social Network
  Analysis: Methods and Applications.} Cambridge: Cambridge University
  Press. }
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
% \seealso{}
\examples{
g <- graph.ring(10)
transitivity(g)
g2 <- erdos.renyi.game(1000, 10/1000)
transitivity(g2)   # this is about 10/1000
}
\keyword{graphs}
