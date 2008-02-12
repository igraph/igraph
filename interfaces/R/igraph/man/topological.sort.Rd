\name{topological.sort}
\alias{topological.sort}
\concept{Topological sort}
\title{Topological sorting of vertices in a graph}
\description{
  A topological sorting of a directed acyclic graph is a linear ordering
  of its nodes where each node comes before all nodes to which it has
  edges. 
}
\usage{
topological.sort(graph, mode=c("out", "all", "in"))
}
\arguments{
  \item{graph}{The input graph, should be directed}
  \item{mode}{Specifies how to use the direction of the edges.
    For \dQuote{\code{out}}, the sorting order ensures that each node comes
    before all nodes to which it has edges, so nodes with no incoming
    edges go first. For \dQuote{\code{in}}, it is quite the opposite: each
    node comes before all nodes from which it receives edges. Nodes 
    with no outgoing edges go first.}
}
\details{
  Every DAG has at least one topological sort, and may have many.
  This function returns a possible topological sort among them. If the
  graph is not acyclic (it has at least one cycle), a partial topological
  sort is returned and a warning is issued.}
\value{
  A numeric vector containing vertex ids in topologically sorted order.
}
%\references{
%}
\author{Tamas Nepusz \email{ntamas@rmki.kfki.hu} and Gabor Csardi
  \email{csardi@rmki.kfki.hu} for the R interface} 
%\seealso{}
\examples{
g <- barabasi.game(100)
topological.sort(g)
}
\keyword{graphs}
