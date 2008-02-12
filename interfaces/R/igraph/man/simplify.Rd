\name{simplify}
\alias{simplify}
\alias{is.simple}
\concept{Simple graph}
\title{Simple graphs}
\description{Simple graphs are graphs which do not contain loop and
  multiple edges.}
\usage{
simplify(graph, remove.multiple = TRUE, remove.loops = TRUE)
is.simple(graph)
}
\arguments{
  \item{graph}{The graph to work on.}
  \item{remove.loops}{Logical, whether the loop edges are to be
    removed.}
  \item{remove.multiple}{Logical, whether the multiple edges are to be
    removed.}
}
\details{
  A loop edge is an edge for which the two endpoints are the same
  vertex. Two edges are multiple edges if they have exactly the same two
  endpoints (for directed graphs order does matter). A graph is simple
  is it does not contain loop edges and multiple edges.

  \code{is.simple} checks whether a graph is simple.

  \code{simplify} removes the loop and/or multiple edges from a graph.
  If both \code{remove.loops} and \code{remove.multiple} are \code{TRUE}
  the function returns a simple graph.
}
\value{A new graph object with the edges deleted. }
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{is.loop}}, \code{\link{is.multiple}} and
  \code{\link{count.multiple}}, \code{\link{delete.edges}},
  \code{\link{delete.vertices}}}
\examples{
g <- graph( c(1,2,1,2,3,3) )
is.simple(g)
is.simple(simplify(g, remove.loops=FALSE))
is.simple(simplify(g, remove.multiple=FALSE))
is.simple(simplify(g))
}
\keyword{graphs}
