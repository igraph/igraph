\name{simplify}
\alias{simplify}
\title{Remove loop and/or multiple edges from a graph}
\description{\code{simplify} removes the multiple and/or loop edges from
a graph.}
\usage{
simplify(graph, remove.loops = TRUE, remove.multiple = TRUE)
}
\arguments{
  \item{graph}{The graph to work on.}
  \item{remove.loops}{Logical, whether the loop edges are to be
    removed.}
  \item{remove.multiple}{Logical, whether the multiple edges are to be
    removed.}
}
\details{If both \code{remove.loops} and \code{remove.multiple} are TRUE
the function returns a simple graph.}
\value{A new graph object with the edges deleted. }
% \references{}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{delete.edges}}, \code{\link{delete.vertices}}}
\examples{
g <- graph( c(1,2,1,2,3,3) )
simplify(g, remove.loops=FALSE)
simplify(g, remove.multiple=FALSE)
}
\keyword{graphs}
