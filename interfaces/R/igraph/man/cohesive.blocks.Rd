\name{cohesive.blocks}
\alias{cohesive.blocks}
\concept{Structurally cohesive blocks}
\title{Calculate Cohesive Blocks}
\description{Calculates cohesive blocks for objects of class \code{igraph}.}
\usage{
cohesive.blocks(graph)
}
\arguments{
    \item{graph}{A graph object of class \code{igraph}. It must de
      undirected and simple. (See \code{\link{is.simple}}.)}
}
\details{
Cohesive blocking is a method of determining hierarchical subsets of
graph vertices based on their structural cohesion (or vertex
connectivity). For a given graph \eqn{G}, a subset of its vertices
\eqn{S\subset V(G)} is said to be maximally \eqn{k}-cohesive if there is
no superset of \eqn{S} with vertex connectivity greater than or equal to
\eqn{k}. Cohesive blocking is a process through which, given a
\eqn{k}-cohesive set of vertices, maximally \eqn{l}-cohesive subsets are
recursively identified with \eqn{l>k}. Thus a hiearchy of vertex subsets
is found, whith the entire graph \eqn{G} at its root. 
}
\value{
  TODO
}
\references{
  J. Moody and D. R. White. Structural cohesion and embeddedness: A
  hierarchical concept of social groups. \emph{American Sociological
    Review}, 68(1):103--127, Feb 2003. 
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com} for the current
  implementation, Peter McMahan \email{peter.mcmahan@gmail.com} wrote
  the first version in R.}
\seealso{\code{\link{graph.cohesion}}}
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
gBlocks$cohesion
}
\keyword{graphs}
