\name{graph.bipartite}
\alias{graph.bipartite}
\concept{Bipartite graph}
\concept{Two-mode network}
\title{Create a bipartite graph}
\description{A bipartite graph has two kinds of vertices and connections
  are only allowed between different kinds.
}
\usage{
graph.bipartite(types, edges, directed=FALSE)
}
\arguments{
  \item{types}{A vector giving the vertex types. It will be coerced into
  boolean. The length of the vector gives the number of vertices in the
  graph.}
  \item{edges}{A vector giving the edges of the graph, the same way as
    for the regular \code{\link{graph}} function. It is checked that the
    edges indeed connect vertices of different kind, accoding to the
    supplied \code{types} vector.}
  \item{directed}{Whether to create a directed graph, boolean
    constant. Note that by default undirected graphs are created, as
    this is more common for bipartite graphs.}
}
\details{
  Bipartite graphs have a \sQuote{\code{type}} vertex attribute in
  igraph, this is boolean and \code{FALSE} for the vertices of the first
  kind and \code{TRUE} for vertices of the second kind.

  \code{graph.bipartite} basically does three things. First it checks
  tha \code{edges} vector against the vertex \code{types}. Then it
  creates a graph using the \code{edges} vector and finally it adds the
  \code{types} vector as a vertex attribute called \code{type}.
}
\value{
  A bipartite igraph graph. In other words, an igraph graph that has a
  vertex attribute \code{type}. 
}
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{graph}} to create one-mode networks}
\examples{
g <- graph.bipartite( rep(0:1,length=10), c(0:9))
print(g, v=TRUE)
}
\keyword{graphs}
