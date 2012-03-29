\name{layout.drl}
\alias{layout.drl}
\alias{igraph.drl.default}
\alias{igraph.drl.coarsen}
\alias{igraph.drl.coarsest}
\alias{igraph.drl.refine}
\alias{igraph.drl.final}
\concept{Graph layout}
\title{The DrL graph layout generator}
\description{
  DrL is a force-directed graph layout toolbox focused on
  real-world large-scale graphs, developed by Shawn Martin and
  colleagues at Sandia National Laboratories.
}
\usage{
layout.drl (graph, use.seed = FALSE, seed = matrix(runif(vcount(graph) * 
    2), ncol = 2), options = igraph.drl.default, weights = E(graph)$weight, 
    fixed = NULL, dim = 2) 
}
\arguments{
  \item{graph}{The input graph, in can be directed or undirected.}
  \item{use.seed}{Logical scalar, whether to use the coordinates given
    in the \code{seed} argument as a starting point.}
  \item{seed}{A matrix with two columns, the starting coordinates for
    the vertices is \code{use.seed} is \code{TRUE}. It is ignored
    otherwise.}
  \item{options}{Options for the layout generator, a named list. See
    details below.}
  \item{weights}{Optional edge weights. Supply \code{NULL} here if you
    want to weight edges equally. By default the \code{weight} edge
    attribute is used if the graph has one.}
  \item{fixed}{Logical vector, it can be used to fix some vertices. All
    vertices for which it is \code{TRUE} are kept at the coordinates
    supplied in the \code{seed} matrix. It is ignored it \code{NULL} or
    if \code{use.seed} is \code{FALSE}. }
  \item{dim}{Either \sQuote{2} or \sQuote{3}, it specifies whether we
    want a two dimensional or a three dimensional layout. Note that
    because of the nature of the DrL algorithm, the three dimensional
    layout takes significantly longer to compute.}
}
\details{
  This function implements the force-directed DrL layout generator.

  The generator has the following parameters:
  \itemize{
    \item{edge.cut}{Edge cutting is done in the late stages of the
      algorithm in order to achieve less dense layouts.  Edges are cut
      if there is a lot of stress on them (a large value in the
      objective function sum). The edge cutting parameter is a value
      between 0 and 1 with 0 representing no edge cutting and 1
      representing maximal edge cutting. } 
    \item{init.iterations}{Number of iterations in the first phase.}
    \item{init.temperature}{Start temperature, first phase.}
    \item{init.attraction}{Attraction, first phase.}
    \item{init.damping.mult}{Damping, first phase.}
    \item{liquid.iterations}{Number of iterations, liquid phase.}
    \item{liquid.temperature}{Start temperature, liquid phase.}
    \item{liquid.attraction}{Attraction, liquid phase.}
    \item{liquid.damping.mult}{Damping, liquid phase.}
    \item{expansion.iterations}{Number of iterations, expansion phase.}
    \item{expansion.temperature}{Start temperature, expansion phase.}
    \item{expansion.attraction}{Attraction, expansion phase.}
    \item{expansion.damping.mult}{Damping, expansion phase.}
    \item{cooldown.iterations}{Number of iterations, cooldown phase.}
    \item{cooldown.temperature}{Start temperature, cooldown phase.}
    \item{cooldown.attraction}{Attraction, cooldown phase.}
    \item{cooldown.damping.mult}{Damping, cooldown phase.}
    \item{crunch.iterations}{Number of iterations, crunch phase.}
    \item{crunch.temperature}{Start temperature, crunch phase.}
    \item{crunch.attraction}{Attraction, crunch phase.}
    \item{crunch.damping.mult}{Damping, crunch phase.}
    \item{simmer.iterations}{Number of iterations, simmer phase.}
    \item{simmer.temperature}{Start temperature, simmer phase.}
    \item{simmer.attraction}{Attraction, simmer phase.}
    \item{simmer.damping.mult}{Damping, simmer phase.}

    There are five pre-defined parameter settings as well, these are
    called \code{igraph.drl.default}, \code{igraph.drl.coarsen},
    \code{igraph.drl.coarsest}, \code{igraph.drl.refine} and
    \code{igraph.drl.final}.
  }
}
\value{
  A numeric matrix with two columns.
}
\references{
  \url{http://www.cs.sandia.gov/~smartin/software.html}
}
\author{Shawn Martin \email{google@for.it} and Gabor Csardi
  \email{csardi@rmki.kfki.hu} for the R/igraph interface and the three
    dimensional version.}
\seealso{\code{\link{layout}} for other layout generators.}
\examples{
g <- as.undirected(ba.game(100, m=1))
l <- layout.drl(g, options=list(simmer.attraction=0))
\dontrun{
plot(g, layout=l, vertex.size=3, vertex.label=NA)
}
}
\keyword{graphs}
