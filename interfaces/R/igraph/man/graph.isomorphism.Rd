\name{graph-isomorphism}
\alias{graph.isoclass}
\alias{graph.isocreate}
\alias{graph.isomorphic}
\alias{graph.isomorphic.vf2}
\title{Graph Isomorphism}
\description{These functions deal with graph isomorphism.}
\usage{
graph.isoclass(graph)
graph.isomorphic(graph1, graph2)
graph.isomorphic.vf2(graph1, graph2)
graph.isocreate(size, number, directed=TRUE)
}
\arguments{
  \item{graph}{A graph object.}
  \item{graph1,graph2}{Graph objects}
  \item{size}{A numeric integer giving the number of vertices in the
    graph to create. Only three or four are suppported right now.}
  \item{number}{The number of the isomorphism class of the graph to be
    created.}
  \item{directed}{Whether to create a directed graph.}
}
\details{
  \code{graph.isoclass} returns the isomorphism class of a graph, a
  non-negative integer number. Graphs (with the same number of vertices)
  having the same isomorphism class are isomorphic and isomorphic graphs
  always have the same isomorphism class.

  \code{graph.isomorphic} decides whether two graphs are isomorphic.

  \code{graph.isomorphic.vf2} decides whethe two graphs are isomorphic,
  it implements the VF2 algorithm, see references.

  \code{graph.isocreate} create a graph from the given isomorphic
  class.
  
  These functions are considered as experimental, as
  \code{graph.isomorphic}, \code{graph.isoclass} and
  \code{graph.isocreate} can handle only graphs with three of four
  vertices.

  It is quite likely that some \code{graph.isomorphic} will be able to
  call \code{graph.isomorphic.vf2} in some later igraph version.
}
\value{
  \code{graph.isoclass} returns a non-negative integer number.

  \code{graph.isomorphic} and \code{graph.isomorphic.vf2}
  return a logical constant.

  \code{graph.isocreate} returns a graph object.
}
\references{
  LP Cordella,  P Foggia, C Sansone, and M Vento:
  An improved algorithm for matching large graphs,
  \emph{Proc. of the 3rd IAPR TC-15 Workshop on Graphbased
    Representations in Pattern Recognition}, 149--159, 2001.
}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{graph.motifs}}}
\examples{
# create some non-isomorphic graphs
g1 <- graph.isocreate(3, 10)
g2 <- graph.isocreate(3, 11)
graph.isoclass(g1)
graph.isomorphic(g1, g2)

# create two isomorphic graphs, by
# permuting the vertices of the first 
g1 <- barabasi.game(30, m=2, directed=FALSE)
el <- get.edgelist(g1)
iso <- sample(vcount(g1))-1
el <- matrix(iso[ el+1 ], nc=2)
g2 <- graph(t(el), directed=FALSE)
graph.isomorphic.vf2(g1, g2)	% should be TRUE
}
\keyword{graphs}
