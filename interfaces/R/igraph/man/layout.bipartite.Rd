\name{layout.bipartite}
\alias{layout.bipartite}
\concept{Graph layout}
\concept{Bipartite graph}
\concept{Two-mode network}
\title{Simple two-row layout for bipartite graphs}
\description{Minimize edge-crossings in a simple two-row (or column)
  layout for bipartite graphs.}
\usage{
layout.bipartite (graph, types = NULL, hgap = 1, vgap = 1, maxiter = 100) 
}
\arguments{
  \item{graph}{The bipartite input graph. It should have a logical
    \sQuote{\code{type}} vertex attribute, or the \code{types} argument
    must be given.}
  \item{types}{A logical vector, the vertex types. If this argument is
    \code{NULL} (the default), then the \sQuote{\code{type}} vertex
    attribute is used.}
  \item{hgap}{Real scalar, the minimum horizontal gap between vertices
    in the same layer.}
  \item{vgap}{Real scalar, the distance between the two layers.}
  \item{maxiter}{Integer scalar, the maximum number of iterations in the
    crossing minimization stage. 100 is a reasonable default; if you
    feel that you have too many edge crossings, increase this.}
}
\details{
  The layout is created by first placing the vertices in two rows,
  according to their types. Then the positions within the rows are
  optimized to minimize edge crossings, using the Sugiyama algorithm
  (see \code{\link{layout.sugiyama}}).
}
\value{A matrix with two columns and as many rows as the number of
  vertices in the input graph.}
%\references{}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{layout.sugiyama}}}
\examples{
# Random bipartite graph
inc <- matrix(sample(0:1, 50, replace=TRUE, prob=c(2,1)), 10, 5)
g <- graph.incidence(inc)
plot(g, layout=layout.bipartite,
     vertex.color=c("green","cyan")[V(g)$type+1])

# Two columns
lay <- layout.bipartite(g)
plot(g, layout=lay[,2:1])
}
\keyword{graphs}
