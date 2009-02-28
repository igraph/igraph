\name{cohesive.blocks}
\alias{cohesive.blocks}
\alias{structurally.cohesive.blocks}
\alias{is.bgraph}
\concept{Structurally cohesive blocks}
\title{Calculate Cohesive Blocks}
\description{Calculates cohesive blocks for objects of class \code{igraph}.}
\usage{
cohesive.blocks(graph, db=NULL,
                useDB=(vcount(graph)>400 && require(RSQLite)),
                cutsetHeuristic=TRUE,
                verbose=igraph.par("verbose"))
is.bgraph(graph)
}
\arguments{
    \item{graph}{A graph object of class \code{igraph}.}
    \item{db}{An optional \code{RSQLite} connection to an existing
      SQLite database (see details). Ignored if \code{NULL}.}
    \item{useDB}{Logical. Whether to use an external SQLite database
      instead of internal R datastructures (see details). By default an
      SQLite database is used if the graph has more than 400 vertices
      and the \code{RSQLite} package is installed.}
    \item{cutsetHeuristic}{Logical scalar. TODO}
    \item{verbose}{Level of console output. Supply \code{TRUE} here to
      follow the progress of the computation.}
}
\details{
Cohesive blocking is a method of determining hierarchical subsets of
graph vertices based on their structural cohesion (or vertex
connectivity). For a given graph \eqn{G}, a subset of its vertices \eqn{S\subset
V(G)} is said to be maximally \eqn{k}-cohesive if there is no superset of
\eqn{S} with vertex connectivity greater than or equal to \eqn{k}. Cohesive
blocking is a process through which, given a \eqn{k}-cohesive set of
vertices, maximally \eqn{l}-cohesive subsets are recursively identified with
\eqn{l>k}. Thus a hiearchy of vertex subsets is found, whith the entire
graph \eqn{G} at its root. 

For certain larger graphs this algorithm can be quite memory-intensive
due to the number of vertex subsets that are examined. In these cases it
is worthwhile to use a database to keep track of this data, specified by
the \code{useDB} argument. If \code{useDB} is \code{TRUE}, then either a
temporary SQLite database is created, or the \code{RSQLite} connection
specified in \code{db} is used. In either case the package
\code{RSQLite} will be required. 

\code{structurally.cohesive.blocks} is an alias to
\code{cohesive.blocks}.

\code{is.bgraph} decides whether its argument is a \code{bgraph}
object.
}
\value{
\code{cohesive.blocks} returns a graph of class \code{c(bgraph,igraph)},
with four (new) graph attributes:
\item{blocks}{A list with one element for each cohesive block found. The
  elements are numeric vectors listing the indices of the nodes within
  that block.} 
\item{block.cohesion}{A numeric vector with length equal to the number
  of cohesive blocks found, listing the cohesion of those blocks.} 
\item{tree}{The hierarchical tree of the cohesive blocks. Each node of
  this graph represents a cohesive block, and directed edges represent
  inclusion as proper subset.} 
\item{data}{A list containing supplementary data from the
  calculation.}

\code{is.bgraph} returns a logical scalar.
}
\references{
  A. Kanevsky. On the number of minimum size separating vertex sets
  in a graph and how to find all of them
  \emph{Proceedings of the first annual ACM-SIAM symposium on Discrete algorithms}
  San Francisco, California, United States. 411--421, 1990.
  
  J. Moody and D. R. White. Structural cohesion and embeddedness: A
  hierarchical concept of social groups. \emph{American Sociological
    Review}, 68(1):103--127, Feb 2003. 
}
\author{Peter McMahan \email{peter.mcmahan@gmail.com}}
\seealso{\code{\link{graph.cohesion}}, \code{\link{plot.bgraph}} for
  plotting graphs together with their block hierarchy,
  \code{\link{write.pajek.bgraph}} for a writing graphs and cohesive
  blocks information to Pajek files. See \code{\link{attributes}} for
  handling graph attributes.} 
\examples{
## Create a graph with an interesting structure:
g <- graph.disjoint.union(graph.full(4), graph.empty(2,directed=FALSE))
g <- add.edges(g,c(3,4,4,5,4,2))
g <- graph.disjoint.union(g,g,g)
g <- add.edges(g,c(0,6,1,7,0,12,4,0,4,1))

## Find cohesive blocks:
gBlocks <- cohesive.blocks(g)

## Examine block membership and cohesion:
gBlocks$blocks
gBlocks$block.cohesion

## Plot the resulting graph with its block hierarchy:
\dontrun{
plot(gBlocks, vertex.size=7, layout=layout.kamada.kawai)
}

## Save the results as Pajek ".net" and ".clu" files:
\dontrun{
write.pajek.bgraph(gBlocks,file="gBlocks")
}

## An example that works better with the "kanevsky" cutset algorithm
\dontrun{
g <- read.graph(file="http://intersci.ss.uci.edu/wiki/Vlado/SanJuanSur.net", format="pajek")
gBlocks <- cohesive.blocks(g, cutsetAlgorithm="kanevsky")
}
}
\keyword{graphs}
