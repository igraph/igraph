\name{bipartite.projection}
\alias{bipartite.projection}
\alias{bipartite.projection.size}
\concept{Bipartite graph}
\concept{Two-mode network}
\title{Project a bipartite graph}
\description{A bipartite graph is projected into two one-mode networks}
\usage{
bipartite.projection.size(graph, types=NULL)
bipartite.projection (graph, types=NULL, multiplicity = TRUE, probe1=-1)
}
\arguments{
  \item{graph}{The input graph. It can be directed, but edge directions
    are ignored during the computation.}
  \item{types}{An optional vertex type vector to use instead of the
    \sQuote{\code{type}} vertex attribute. You must supply this argument
    if the graph has no \sQuote{\code{type}} vertex attribute.}
  \item{multiplicity}{If \code{TRUE}, then igraph keeps the multiplicity
    of the edges as an edge attribute. E.g. if there is an A-C-B and
    also an A-D-B triple in the bipartite graph (but no more X, such
    that A-X-B is also in the graph), then the multiplicity of the A-B
    edge in the projection will be 2.}
  \item{probe1}{This argument can be used to specify the order of the
    projections in the resulting list. If given and non-negative, then
    it is considered as a vertex id; the projection containing the
    first one in the result list.}
}
\details{
  Bipartite graphs have a \sQuote{\code{type}} vertex attribute in
  igraph, this is boolean and \code{FALSE} for the vertices of the first
  kind and \code{TRUE} for vertices of the second kind.

  \code{bipartite.projection.size} calculates the number of vertices and
  edges in the two projections of the bipartite graphs, without
  calculating the projections themselves. This is useful to check how
  much memory the projections would need if you have a large bipartite
  graph.

  \code{bipartite.projections} calculates the actual projections.
  You can use the \code{probe1} argument to specify the order of the
  projections in the result. By default vertex type \code{FALSE} is the
  first and \code{TRUE} is the second.

  \code{bipartite.projections} keeps vertex attributes.
}
\value{
  A list of two undirected graphs. See details above.
}
% \references{}
\author{ Gabor Csardi \email{csardi@rmki.kfki.hu} }
% \seealso{}
\examples{
## Projection of a full bipartite graph is a full graph
g <- graph.full.bipartite(10,5)
proj <- bipartite.projection(g)
graph.isomorphic(proj[[1]], graph.full(10))
graph.isomorphic(proj[[2]], graph.full(5))

## The projection keeps the vertex attributes
M <- matrix(0, nr=5, nc=3)
rownames(M) <- c("Alice", "Bob", "Cecil", "Dan", "Ethel")
colnames(M) <- c("Party", "Skiing", "Badminton")
M[] <- sample(0:1, length(M), replace=TRUE)
M
g2 <- graph.incidence(M)
g2$name <- "Event network"
proj2 <- bipartite.projection(g2)
print(proj2[[1]], g=TRUE, e=TRUE)
print(proj2[[2]], g=TRUE, e=TRUE)
}
\keyword{graphs}
