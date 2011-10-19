\name{plot.bgraph}
\alias{plot.bgraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot graphs and their cohesive block hierarchy}
\description{
  This function plots \code{bgraph} objects as output by
  \code{cohesive.blocks}. It produces a two-panel plot with the graph
  itself on the left and a tree representing the block hierarchy on the
  right. 
}
\usage{
\method{plot}{bgraph}(x, mc = NULL, vertex.size = 3, colpal = NULL, emph = NULL, ...)
}
\arguments{
  \item{x}{The \code{bgraph} object to be plotted}
  \item{mc}{A numeric vector listing the vertex connectivity of the
    maximally cohesive subgraph of each vertex. Automatically calculated
    if \code{NULL} (leaving default is usually preferable).} 
  \item{vertex.size}{The size of the vertices in the plot. Applies only
    to the graph, and not to the block-hierarchy tree.} 
  \item{colpal}{The color palette to use to distinguish block
    cohesion. Defaults to rainbow spectrum.} 
  \item{emph}{A numeric vector of blocks to emphasize. Useful for
    distinguishing specific blocks when it is unclear which
    higher-cohesion vertices belong to which block. (see details)} 
  \item{\dots}{Other arguments to be passed on to \code{plot.igraph} for
    the calculation of the graph (but not the hierarchy).}
}
\details{
Two plots are used to represent the cohesive blocks in a graph
visually. The first is a standard plot with vertices colored according
to their maximally-cohesive containing block. The second is a tree
representing the hierarchical structure of the blocks, with edges
representing a strict superset relationship.

The \code{emph} argument should be a numeric vector corresponding to the
indices of blocks in \code{x$blocks} and \code{x$block.cohesion}
(1-based indexing). The vertices of the specified blocks are emphasized
by enlarging them and using a white border.

The intended usage of this function is the quick plotting of a graph
together with its block structure. If you need more flexibility then
please plot the graph and the hierarchy (the \code{tree} graph
attribute) separately by using \code{\link{plot.igraph}}.
}
\author{Peter McMahan \email{peter.mcmahan@gmail.com}}
\seealso{\code{\link{cohesive.blocks}} for the cohesive blocks
  computation, \code{\link{graph.cohesion}}, \code{\link{plot.igraph}}
  and \code{\link{igraph.plotting}} for regular igraph plotting,
  \code{\link{write.pajek.bgraph}}.}
\examples{
## Create a graph with an interesting structure:
g <- graph.disjoint.union(graph.full(4),graph.empty(2,directed=FALSE))
g <- add.edges(g,c(3,4,4,5,4,2))
g <- graph.disjoint.union(g,g,g)
g <- add.edges(g,c(0,6,1,7,0,12,4,0,4,1))

## Find cohesive blocks:
gBlocks <- cohesive.blocks(g)

## Plot:
\dontrun{
plot.bgraph(gBlocks,layout=layout.kamada.kawai)
}

## There are two two-cohesive blocks. To differentiate the block 
## that contains both the three- and four-cohesive sub-blocks use:
\dontrun{
plot(gBlocks,emph=3,layout=layout.kamada.kawai)
}
}
\keyword{graphs}
