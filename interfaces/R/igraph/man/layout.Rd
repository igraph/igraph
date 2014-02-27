\name{layout}
\alias{layout}
\alias{layout.auto}
\alias{layout.random}
\alias{layout.circle}
\alias{layout.sphere}
\alias{layout.spring}
\alias{layout.reingold.tilford}
\alias{layout.lgl}
\alias{layout.svd}
\alias{layout.graphopt}
\alias{layout.norm}
\concept{Graph layout}
\title{Generate coordinates for plotting graphs}
\description{Some simple and not so simple functions determining the
  placement of the vertices for drawing a graph.}
\usage{
layout.auto(graph, dim=2, \dots)
layout.random(graph, params, dim=2)
layout.circle(graph, params)
layout.sphere(graph, params)
layout.spring(graph, \dots, params)
layout.reingold.tilford(graph, \dots, params)
layout.lgl(graph, \dots, params)
layout.graphopt(graph, ..., params=list())
layout.svd(graph, d=shortest.paths(graph), ...)
layout.norm(layout, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL,
      zmin = NULL, zmax = NULL)
}
\arguments{
  \item{graph}{The graph to place.}
  \item{params}{The list of function dependent parameters.}
  \item{dim}{Numeric constant, either 2 or 3. Some functions are able
    to generate 2d and 3d layouts as well, supply this argument to
    change the default behavior.}
  \item{\dots}{Function dependent parameters, this is an alternative
    notation to the \code{params} argument. For \code{layout.auto} these
    extra parameters are simply passed to the real layout function, if
    one is called.}
  \item{d}{The matrix used for singular value decomposition. By default it
    is the distance matrix of the graph.}
  \item{layout}{A matrix with two or three columns, the layout to
    normalize.}
  \item{xmin,xmax}{The limits for the first coordinate, if one of them
    or both are \code{NULL} then no normalization is performed along
    this direction.}
  \item{ymin,ymax}{The limits for the second coordinate, if one of them
    or both are \code{NULL} then no normalization is performed along
    this direction.}
  \item{zmin,zmax}{The limits for the third coordinate, if one of them
    or both are \code{NULL} then no normalization is performed along
    this direction.}  
}
\details{
  These functions calculate the coordinates of the vertices for a graph
  usually based on some optimality criterion.

  \code{layout.auto} tries to choose an appropriate layout function for
  the supplied graph, and uses that to generate the layout. The current
  implementations works like this:
  \enumerate{
    \item If the graph has a graph attribute called \sQuote{layout},
    then this is used. If this attribute is an R function, then it is
    called, with the graph and any other extra arguments.
    \item Otherwise, if the graph has vertex attributes called
    \sQuote{x} and \sQuote{y}, then these are used as coordinates. If
    the graph has an additional \sQuote{z} vertex attribute, that is
    also used.
    \item Otherwise, if the graph is connected and has less than 100
    vertices, the Kamada-Kawai layout is used, by calling
    \code{layout.kamada.kawai}.
    \item Otherwise, if the graph has less than 1000 vertices, then the
    Fruchterman-Reingold layout is used, by calling
    \code{layout.fruchterman.reingold}.
    \item Otherwise the DrL layout is used, \code{layout.drl} is called.
  }
  
  \code{layout.random} simply places the vertices randomly on a
  square. It has no parameters.

  \code{layout.circle} places the vertices on a unit circle
  equidistantly. It has no paramaters.

  \code{layout.sphere} places the vertices (approximately) uniformly on
  the surface of a sphere, this is thus a 3d layout. It is not clear
  however what \dQuote{uniformly on a sphere} means.
  
  \code{layout.spring} is a spring embedder algorithm.
  Parameters and default values:
  \describe{
    \item{mass}{The vertex mass (in \sQuote{quasi-kilograms}). (Defaults to 0.1.)}
    \item{equil}{The equilibrium spring extension (in
      \sQuote{quasi-meters}). (Defaults to 1.)}
    \item{k}{The spring coefficient (in \sQuote{quasi-Newtons per
      quasi-meter}). (Defaults to 0.001.)}
    \item{repeqdis}{The point at which repulsion (if employed) balances
      out the spring extension force (in \sQuote{quasi-meters}). (Defaults to
      0.1.)}
    \item{kfr}{The base coefficient of kinetic friction (in
      \sQuote{quasi-Newton quasi-kilograms}). (Defaults to 0.01.)}
    \item{repulse}{Should repulsion be used? (Defaults to FALSE.)}
  }
  This function was ported from the SNA package.

  \code{layout.reingold.tilford} generates a tree-like layout, so it is
  mainly for trees. Parameters and default values:
  \describe{
    \item{root}{The id of the root vertex, defaults to 1.}
    \item{circular}{Logical scalar, whether to plot the tree in a
      circular fashion, defaults to \code{FALSE}.}
    \item{flip.y}{Logical scalar, whether to flip the \sQuote{y}
      coordinates. The default is flipping because that puts the root
      vertex on the top.}
  }
  
  \code{layout.lgl} is for large connected graphs, it is similar to the
  layout generator of the Large Graph Layout software
  (\url{http://lgl.sourceforge.net/}). Parameters and
  default values:
  \describe{
    \item{maxiter}{The maximum number of iterations to perform (150).}
    \item{maxdelta}{The maximum change for a vertex during an iteration
      (the number of vertices).}
    \item{area}{The area of the surface on which the vertices are placed
      (square of the number of vertices). }
    \item{coolexp}{The cooling exponent of the simulated annealing (1.5).}
    \item{repulserad}{Cancellation radius for the repulsion (the
      \code{area} times the number of vertices).}
    \item{cellsize}{The size of the cells for the grid. When calculating
      the repulsion forces between vertices only vertices in the same or
      neighboring grid cells are taken into account (the fourth root of
      the number of \code{area}.}
    \item{root}{The id of the vertex to place at the middle of the
      layout. The default value is -1 which means that a random vertex
      is selected.}
  }

  \code{layout.graphopt} is a port of the graphopt layout algorithm by
  Michael Schmuhl. graphopt version 0.4.1 was rewritten in C and the
  support for layers was removed (might be added later) and a code was a
  bit reorganized to avoid some unneccessary steps is the node charge
  (see below) is zero.

  graphopt uses physical analogies for defining attracting and repelling 
  forces among the vertices and then the physical system is simulated 
  until it reaches an equilibrium. (There is no simulated annealing or 
  anything like that, so a stable fixed point is not guaranteed.)
  
  See also \url{http://www.schmuhl.org/graphopt/} for the original
  graphopt.

  Parameters and default values:
  \describe{
    \item{niter}{Integer scalar, the number of iterations to perform.
      Should be a couple of hundred in general. If you have a large graph 
      then you might want to only do a few iterations and then check the 
      result. If it is not good enough you can feed it in again in 
      the \code{start} argument. The default value is 500.
    }
    \item{charge}{The charge of the vertices, used to calculate electric
      repulsion. The default is 0.001.}
    \item{mass}{The mass of the vertices, used for the spring
      forces. The default is 30.}
    \item{spring.length}{The length of the springs, an integer
      number. The default value is zero.}
    \item{spring.constant}{The spring constant, the default value is one.}
    \item{max.sa.movement}{Real constant, it gives the maximum amount of
      movement allowed in a single step along a single axis. The default
      value is 5.}
    \item{start}{If given, then it should be a matrix with two columns
      and one line for each vertex. This matrix will be used as starting
      positions for the algorithm. If not given, then a random starting
      matrix is used.}
  }
  
  \code{layout.svd} is a currently \emph{experimental} layout function based on
  singular value decomposition. It does not have the usual \code{params}
  argument, but take a single argument, the distance matrix of the
  graph. This function generates the layout separately for each graph
  component and then merges them via \code{\link{layout.merge}}.

  \code{layout.norm} normalizes a layout, it linearly transforms each
  coordinate separately to fit into the given limits.

  \code{layout.drl} is another force-driven layout generator, it is
  suitable for quite large graphs. See \code{\link{layout.drl}} for
  details.
}
\value{All these functions return a numeric matrix with at least two
  columns and the same number of lines as the number of vertices.}
\references{
  Reingold, E and Tilford, J (1981). Tidier drawing of trees.
  \emph{IEEE Trans. on Softw. Eng.}, SE-7(2):223--228.
}
\author{Gabor Csardi \email{csardi.gabor@gmail.com}}
\seealso{\code{\link{layout.drl}}, \code{\link{plot.igraph}}, \code{\link{tkplot}}}
\examples{
g <- graph.ring(10)
layout.random(g)
}
\keyword{graphs}
