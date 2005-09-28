\name{layout}
\alias{layout}
\alias{layout.random}
\alias{layout.circle}
\alias{layout.fruchterman.reingold}
\alias{layout.kamada.kawai}
\alias{layout.spring}
\title{Generate coordinates for plotting graphs}
\description{Some simple and not so simple functions determining the
  placement of the vertices for drawing a graph.}
\usage{
layout.random(graph, params)
layout.circle(graph, params)
layout.fruchterman.reingold(graph, \dots, params)
layout.kamada.kawai(graph, \dots, params)
layout.spring(graph, \dots, params)
}
\arguments{
  \item{graph}{The graph to place.}
  \item{params}{The list of function dependent parameters.}
  \item{\dots}{Function dependent parameters, this is an alternative
    notation to the \code{params} argument.}
}
\details{
  These functions calculate the coordinates of the vertices for a graph
  usually based on some optimality criterion.

  \code{layout.random} simply places the vertices randomly on a
  square. It has no parameters.

  \code{layout.circle} places the vertices on a unit circle
  equidistantly. It has no paramaters.

  \code{layout.fruchterman.reingold} uses a force-based algorithm
  proposed by Fruchterman and Reingold, see references. Parameters and
  their default values:
  \itemize{
    \item{niter}{Numeric, the number of iterations to perform (100).}
    \item{coolexp}{Numeric, the cooling exponent for the simulated
      annealing (1.5).}
    \item{frame}{Character, the type of the frame to place the
      vertices. It is either \code{rectangle} or \code{circle}
      (rectangle). }
    \item{initial}{Numeric matrix, the layout to start with (random).}
    \item{temp}{The initial temperature for the simulated annealing
      (0.1).} 
  }

  \code{layout.kamada.kawai} is other force based algorithm. This was
  ported from the SNA package. Parameters and default values:
  \itemize{
    \item{niter}{Number of iterations to perform (1000).}
    \item{sigma}{Sets the base standard deviation of position change
      proposals (vcount(graph)/4).}
    \item{initemp}{The initial temperature (10).}
    \item{coolexp}{The cooling exponent (0.99).}
    \item{kkconst}{Sets the Kamada-Kawai vertex attraction
      constant (vcount(graph)**2).}
  }

  \code{layout.spring} is a spring embedder algorithm. This function was
  ported from the SNA package. Parameters and default values:
  \itemize{
    \item{mass}{The vertex mass (in ``quasi-kilograms''). (Defaults to 0.1.)}
    \item{equil}{The equilibrium spring extension (in
      ``quasi-meters''). (Defaults to 1.)}
    \item{k}{The spring coefficient (in ``quasi-Newtons per
      quasi-meter''). (Defaults to 0.001.)}
    \item{repeqdis}{The point at which repulsion (if employed) balances
      out the spring extension force (in ``quasi-meters''). (Defaults to
      0.1.)}
    \item{kfr}{The base coefficient of kinetic friction (in
      ``quasi-Newton quasi-kilograms''). (Defaults to 0.01.)}
    \item{repulse}{Should repulsion be used? (Defaults to FALSE.)}
  }    
}
\value{All these functions return a numeric matrix with at least two
  columns and the same number of lines as the number of vertices.}
\references{
  Fruchterman, T.M.J. and Reingold, E.M. (1991). ``Graph Drawing by
  Force-directed Placement.'' Software - Practice and Experience,
  21(11):1129-1164. 

  Kamada, T. and Kawai, S. (1989). ``An Algorithm for Drawing General
  Undirected Graphs.'' Information Processing Letters, 31(1):7-15. 
}
\author{Gabor Csardi \email{csardi@rmki.kfki.hu}}
\seealso{\code{\link{plot.igraph}}, \code{\link{tkplot}}}
\examples{
g <- graph.ring(10)
layout.random(g)
layout.kamada.kawai(g)
}
\keyword{graphs}
