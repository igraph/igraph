\name{graph.density}
\alias{graph.density}
\concept{Graph density}
\title{Graph density}
\description{The density of a graph is the ratio of the number of edges
  and the number of possible edges.}
\usage{
graph.density(graph, loops=FALSE)
}
\arguments{
  \item{graph}{The input graph.}
  \item{loops}{Logical constant, whether to allow loop edges in the
    graph. If this is TRUE then self loops are considered to be possible. If
    this is FALSE then we assume that the graph does not contain any
    loop edges and that loop edges are not meaningful.}
}
\details{
  Note that this function may return strange results for graph with
  multiple edges, density is ill-defined for graphs with multiple
  edges. 
}
\value{
  A real constant. This function returns \code{NaN} (=0.0/0.0) for an
  empty graph with zero vertices.
}
\references{
  Wasserman, S., and Faust, K.  (1994).  Social Network Analysis:
  Methods and Applications.  Cambridge: Cambridge University Press.
}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{vcount}}, \code{\link{ecount}},
  \code{\link{simplify}} to get rid of the multiple and/or loop edges.}
\examples{
g1 <- graph.empty(n=10)
g2 <- graph.full(n=10)
g3 <- erdos.renyi.game(n=10, 0.4)

# loop edges
g <- graph( c(0,1, 1,1, 1,2) )
graph.density(g, loops=FALSE)              # this is wrong!!!
graph.density(g, loops=TRUE)               # this is right!!!
graph.density(simplify(g), loops=FALSE)    # this is also right, but different
}
\keyword{graphs}
