\name{graph.kcores}
\alias{graph.kcores}
\title{K-core decomposition of graphs}
\description{The k-core of graph is a maximal subgraph in which each
  vertex has at least degree k.}
\usage{
graph.kcores(graph, mode="all")
}
\arguments{
  \item{graph}{The input graph, it can be directed or undirected}
  \item{mode}{The type of the core in directed graphs. Character
    constant, possible values: \code{in}: in-cores are computed,
    \code{out}: out-cores are computed, \code{all}: the corresponding
    undirected graph is considered. This argument is ignored for
    undirected graphs.}
}
\details{
  The k-core of a graph is the maximal subgraph in which every vertex
  has at least degree k. The cores of a graph form layers: the (k+1)-core
  is always a subgraph of the k-core.

  This function calculates the order for each vertex: the
  maximum order of a core in which it is included.
}
\value{
  Numeric vector of integer numbers giving the order of each vertex.
}
\references{Vladimir Batagelj, Matjaz Zaversnik: An O(m)
  Algorithm for Cores Decomposition of Networks, 2002
  
  Seidman S. B. (1983) Network structure and minimum degree, \emph{Social
    Networks}, 5, 269--287.

}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{degree}}}
\examples{
g <- graph.ring(10)
g <- add.edges(g, c(0,1, 1,2, 0,2))
graph.kcores(g) 		% small core triangle in a ring
}
\keyword{graphs}