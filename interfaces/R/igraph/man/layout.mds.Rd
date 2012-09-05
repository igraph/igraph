\name{layout.mds}
\alias{layout.mds}
\concept{Graph layout}
\title{Graph layout by multidimensional scaling}
\description{
  Multidimensional scaling of some distance matrix defined on the
  vertices of a graph.
}
\usage{
layout.mds(graph, dist=NULL, dim=2, options=igraph.arpack.default,
           returnParams=FALSE)
}
\arguments{
  \item{graph}{The input graph.}
  \item{dist}{The distance matrix for the multidimensional scaling.
    If \code{NULL} (the default), then the unweighted shortest path
    matrix is used.}
  \item{dim}{\code{layout.mds} supports dimensions up to the number of
    nodes minus one, but only if the graph is connected; for unconnected
    graphs, the only possible values is 2. This is because
    \code{layout.merge} only works in 2D.}
  \item{options}{A named list, to override some ARPACK options. See
    \code{\link{arpack}} for details.}
  \item{returnParams}{Logical scalar. Whether to return ARPACK
    parameters from the fitting, or just the layout coordinates in a
    matrix. By default only the matrix is returned.}
}
\details{
  \code{layout.mds} uses metric multidimensional scaling for generating
  the coordinates. Multidimensional scaling aims to place points from a
  higher dimensional space in a (typically) 2 dimensional plane, so that
  the distance between the points are kept as much as this is possible.

  By default igraph uses the shortest path matrix as the distances
  between the nodes, but the user can override this via the \code{dist}
  argument.

  This function generates the layout separately for each graph component
  and then merges them via \code{\link{layout.merge}}.
}
\value{
  Depending on the value of the \code{returnParams}, this function
  returns either a numeric matrix with \code{dim} columns, or a list
  with two entries:
  \item{res}{The layout matrix.}
  \item{options}{The parameters from the ARPACK run. See
  \code{\link{arpack}} for details.}
}
\references{
  Cox, T. F. and Cox, M. A. A. (2001) \emph{Multidimensional Scaling}.
  Second edition. Chapman and Hall.
}
\author{Tamas Nepusz \email{ntamas@gmail.com} and
  Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{layout}}, \code{\link{plot.igraph}}}
\examples{
g <- erdos.renyi.game(100, 2/100)
l <- layout.mds(g)
plot(g, layout=l, vertex.label=NA, vertex.size=3)
}
\keyword{graphs}
