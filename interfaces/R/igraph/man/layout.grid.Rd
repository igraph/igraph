\name{layout.grid}
\alias{layout.grid}
\alias{layout.grid.3d}
\concept{Graph layout}
\title{Simple grid layout}
\description{This layout places vertices on a rectangulat grid, in two
  or three dimensions.}
\usage{
layout.grid (graph, width = 0)
layout.grid.3d (graph, width = 0, height = 0)
}
\arguments{
  \item{graph}{The input graph.}
  \item{width}{The number of vertices in a single row of the grid. If
    this is zero or negative for \code{layout.grid}, then the width of
    the grid will be the square root of the number of vertices in the
    graph, rounded up to the next integer. Similarly, it will be the
    cube root for \code{layout.grid.3d}.}
  \item{height}{The number of vertices in a single column of the
    grid, for three dimensional layouts. If this is zero or negative,
    then it is determinted automatically.}
}
\details{
  These functions place the vertices on a simple rectangular grid, one
  after the other. If you want to change the order of the vertices, then
  see the \code{\link{permute.vertices}} function.
}
\value{
  A two-column matrix for \code{layout.grid}, a three-column matrix for
  \code{layout.grid.3d}.
}
% \references{}
\author{ Tamas Nepusz \email{ntamas@gmail.com} }
\seealso{ \code{\link{layout}} for other layout generators }
\examples{
g <- graph.lattice( c(3,3) )
layout.grid(g)

g2 <- graph.lattice( c(3,3,3) )
layout.grid.3d(g2)

\dontrun{
plot(g, layout=layout.grid)
rglplot(g, layout=layout.grid.3d)
}
}
\keyword{graphs}
