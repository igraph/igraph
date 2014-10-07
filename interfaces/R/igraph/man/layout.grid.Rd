\name{layout_on_grid}
\alias{layout_on_grid}
\alias{layout.grid}
\alias{layout.grid.3d}
\concept{Graph layout}
\title{Simple grid layout}
\description{This layout places vertices on a rectangulat grid, in two
  or three dimensions.}
\usage{
layout_on_grid (graph, width = 0, height = 0, dim = 2)
}
\arguments{
  \item{graph}{The input graph.}
  \item{width}{The number of vertices in a single row of the grid. If
    this is zero or negative, then for 2d layouts the width of the grid
    will be the square root of the number of vertices in the graph,
    rounded up to the next integer. Similarly, it will be the 
    cube root for 3d layouts.}
  \item{height}{The number of vertices in a single column of the
    grid, for three dimensional layouts. If this is zero or negative,
    then it is determinted automatically.}
  \item{dim}{Two or three. Whether to make 2d or a 3d layout.}
}
\details{
  The function places the vertices on a simple rectangular grid, one
  after the other. If you want to change the order of the vertices, then
  see the \code{\link{permute}} function.
}
\value{
  A two-column or three-column matrix.
}
% \references{}
\author{ Tamas Nepusz \email{ntamas@gmail.com} }
\seealso{ \code{\link{layout}} for other layout generators }
\examples{
g <- lattice( c(3,3) )
layout_on_grid(g)

g2 <- lattice( c(3,3,3) )
layout_on_grid(g2, dim = 3)

\dontrun{
plot(g, layout=layout_on_grid)
rglplot(g, layout=layout_on_grid(g, dim = 3))
}
}
\keyword{graphs}
