\name{layout.graphopt}
\alias{layout.graphopt}
\concept{Graph layout}
\title{The graphopt layout algorithm}
\description{
  A force-directed layout algorithm, that scales relatively well to
  large graphs.
}
\usage{
layout.graphopt (graph, start=NULL, niter=500, charge=0.001,
                 mass=30, spring.length=0, spring.constant=1,
                 max.sa.movement=5)
}
\arguments{
  \item{graph}{The input graph.}
  \item{start}{If given, then it should be a matrix with two columns
    and one line for each vertex. This matrix will be used as starting
    positions for the algorithm. If not given, then a random starting
    matrix is used.}
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
}
\details{
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
}
\value{
  A numeric matrix with two columns, and a row for each vertex.
}
% \references{}
\author{Michael Schmuhl for the original graphopt code,
  rewritten and wrapped by Gabor Csardi \email{csardi.gabor@gmail.com}.
}
% \seealso{}
% \examples{}
\keyword{graphs}

