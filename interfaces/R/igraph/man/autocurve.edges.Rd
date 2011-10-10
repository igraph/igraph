\name{autocurve.edges}
\alias{autocurve.edges}
\title{Optimal edge curvature when plotting graphs}
\description{
  If graphs have multiple edges, then drawing them as straight lines
  does not show them when plotting the graphs; they will be on top
  of each other. One solution is to bend the edges, with diffenent
  curvature, so that all of them are visible.}
\usage{
autocurve.edges (graph, start = 0.5)
}
\arguments{
  \item{graph}{The input graph.}
  \item{start}{The curvature at the two extreme edges. All edges will
    have a curvature between \code{-start} and \code{start}, spaced
    equally.}
}
\details{
  \code{autocurve.edges} calculates the optimal \code{edge.curved}
  vector for plotting a graph with multiple edges, so that all edges
  are visible.
}
\value{
  A numeric vector, its length is the number of edges in the graph.
}
% \references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{igraph.plotting}} for all plotting parameters,
  \code{\link{plot.igraph}}, \code{\link{tkplot}} and
  \code{\link{rglplot}} for plotting functions.
}
\examples{
g <- graph( c(0,1,1,0,1,2,1,3,1,3,1,3,
              2,3,2,3,2,3,2,3,0,1)+1 )

autocurve.edges(g)

\dontrun{
set.seed(42)
plot(g)
}
}
\keyword{graphs}
