\name{edge.betweenness.community}
\alias{edge.betweenness.community}
\alias{edge.betweenness.community.merges}
\concept{Edge betweenness}
\concept{Community structure}
\title{Community structure detection based on edge betweenness}
\description{Many networks consist of modules which are densely
  connected themselves but sparsely connected to other modules.}
\usage{
edge.betweenness.community (graph, directed = TRUE,
    edge.betweenness = TRUE, merges = TRUE, bridges = TRUE,
    labels = TRUE)
edge.betweenness.community.merges (graph, edges)
}
\arguments{
  \item{graph}{The graph to analyze.}
  \item{directed}{Logical constant, whether to calculate directed edge
    betweenness for directed graphs. It is ignored for undirected
    graphs.}
  \item{edge.betweenness}{Logical constant, whether to return the edge
    betweenness of the edges at the time of their removal.}
  \item{merges}{Logical constant, whether to return the merge matrix
    representing the hierarchical community structure of the network.
    This argument is called \code{merges}, even if the community
    structure algorithm itself is divisive and not agglomerative: it
    builds the tree from top to bottom. There is one line for each merge
    (ie. split) in matrix, the first line is the first merge (last
    split). The communities are identified by integer number starting from
    zero. Community ids smaller than \sQuote{\code{N}}, the number of
    vertices in the graph, belong to singleton communities, ie. individual
    vertices. Before the first merge we have \sQuote{\code{N}} communities
    numbered from zero to \sQuote{\code{N-1}}. The first merge, the
    first line of the matrix creates community \sQuote{\code{N}}, the
    second merge creates community \sQuote{\code{N+1}}, etc.
  }
  \item{bridges}{Logical constant, whether to return a list the edge
    removals which actually splitted a component of the graph.}
  \item{labels}{Logical constant, whether to contain the labels of the
    vertices in the result. More precisely, if the graph has a vertex
    attribute valled \sQuote{name}, it will be part of the result
    object.}
  \item{edges}{Numeric vector, the ids of the edges to be removed from a
    graph, all edges should be present in the vector, their order
    specifies the order of removal.}
}
\details{
  The edge betweenness score of an edge measures the number of shortest
  paths through it, see \code{\link{edge.betweenness}} for details. The
  idea of the edge betweenness based community structure detection is
  that it is likely that edges connecting separate modules have high
  edge betweenness as all the shortest paths from one module to another
  must traverse through them. So if we gradually remove the edge with
  the highest edge betweenness score we will get a hierarchical map, a
  rooted tree, called a dendrogram of the graph. The leafs of the tree
  are the individual vertices and the root of the tree represents the
  whole graph. 

  \code{edge.betweenness.community} performs this algorithm by
  calculating the edge betweenness of the graph, removing the edge with
  the highest edge betweenness score, then recalculating edge
  betweenness of the edges and again removing the one with the highest
  score, etc.

  \code{edge.betweeness.community} returns various information collected
  throught the run of the algorithm. See the return value down here.

  \code{edge.betweenness.community.merges} gets a list of edges and by
  gradually removes them from the graph it creates a merge matrix
  similar to the one returned by \code{edge.betweenness.community}.
}
\value{
  A named list is returned by \code{edge.betweenness.community},
  with the following components:
  \item{removed.edges}{Numeric vector, the edges of the graph, in the
    order of their removal.}
  \item{edge.betweenness}{Numeric vector, the edge betweenness value
    of the removed edges, the order is the same as in
    \code{removed.edges}.}
  \item{merges}{Matrix containing the merges (ie. divisions) the
    algorithm performed, see the \code{merges} argument for the
    format.}
  \item{bridges}{Numeric vector, the steps (ie. edge removals) which
    resulted a split of a component in the graph.} 
  \item{labels}{The \code{name} argument of the vertices.}

  Note that some components may be missing or \code{NULL} if you do not
  request them, see the parameters.

  A numeric matrix is returned by
  \code{edge.betweenness.community.merges}.
  The matrix has two column and its format is the same as the
  \code{merges} slot of the result of
  \code{edge.betweenness.community}.
}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{ \code{\link{edge.betweenness}} for the definition and
  calculation of the edge betweenness,
  \code{\link{walktrap.community}}, \code{\link{fastgreedy.community}},
  \code{\link{leading.eigenvector.community}} for other community
  detection methods. \code{as.dendrogram} in package \code{stats} for creating an
  R dendrogram object from the result of the clustering. See
  \code{\link{community.to.membership}} to create the actual communities
  after a number of edges removed from the network.
}
\references{M Newman and M Girvan: Finding and evaluating community
  structure in networks, \emph{Physical Review E} 69, 026113 (2004)
}
\examples{
g <- barabasi.game(100,m=2)
eb <- edge.betweenness.community(g)

g <- graph.full(10) \%du\% graph.full(10)
g <- add.edges(g, c(0,10))
eb <- edge.betweenness.community(g)
E(g) [ eb$removed.edges[1] ] 		% should be the edge just added
}
\keyword{graphs}
