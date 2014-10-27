
#' The DrL graph layout generator
#' 
#' DrL is a force-directed graph layout toolbox focused on real-world
#' large-scale graphs, developed by Shawn Martin and colleagues at Sandia
#' National Laboratories.
#' 
#' This function implements the force-directed DrL layout generator.
#' 
#' The generator has the following parameters: \describe{ \item{edge.cut}{Edge
#' cutting is done in the late stages of the algorithm in order to achieve less
#' dense layouts.  Edges are cut if there is a lot of stress on them (a large
#' value in the objective function sum). The edge cutting parameter is a value
#' between 0 and 1 with 0 representing no edge cutting and 1 representing
#' maximal edge cutting. } \item{init.iterations}{Number of iterations in the
#' first phase.} \item{init.temperature}{Start temperature, first phase.}
#' \item{init.attraction}{Attraction, first phase.}
#' \item{init.damping.mult}{Damping, first phase.}
#' \item{liquid.iterations}{Number of iterations, liquid phase.}
#' \item{liquid.temperature}{Start temperature, liquid phase.}
#' \item{liquid.attraction}{Attraction, liquid phase.}
#' \item{liquid.damping.mult}{Damping, liquid phase.}
#' \item{expansion.iterations}{Number of iterations, expansion phase.}
#' \item{expansion.temperature}{Start temperature, expansion phase.}
#' \item{expansion.attraction}{Attraction, expansion phase.}
#' \item{expansion.damping.mult}{Damping, expansion phase.}
#' \item{cooldown.iterations}{Number of iterations, cooldown phase.}
#' \item{cooldown.temperature}{Start temperature, cooldown phase.}
#' \item{cooldown.attraction}{Attraction, cooldown phase.}
#' \item{cooldown.damping.mult}{Damping, cooldown phase.}
#' \item{crunch.iterations}{Number of iterations, crunch phase.}
#' \item{crunch.temperature}{Start temperature, crunch phase.}
#' \item{crunch.attraction}{Attraction, crunch phase.}
#' \item{crunch.damping.mult}{Damping, crunch phase.}
#' \item{simmer.iterations}{Number of iterations, simmer phase.}
#' \item{simmer.temperature}{Start temperature, simmer phase.}
#' \item{simmer.attraction}{Attraction, simmer phase.}
#' \item{simmer.damping.mult}{Damping, simmer phase.}
#' 
#' There are five pre-defined parameter settings as well, these are called
#' \code{drl_defaults$default}, \code{drl_defaults$coarsen},
#' \code{drl_defaults$coarsest}, \code{drl_defaults$refine} and
#' \code{drl_defaults$final}.  }
#' 
#' @aliases layout.drl drl_defaults igraph.drl.coarsen
#'  igraph.drl.coarsest igraph.drl.default igraph.drl.final
#'  igraph.drl.refine
#' @param graph The input graph, in can be directed or undirected.
#' @param use.seed Logical scalar, whether to use the coordinates given in the
#' \code{seed} argument as a starting point.
#' @param seed A matrix with two columns, the starting coordinates for the
#' vertices is \code{use.seed} is \code{TRUE}. It is ignored otherwise.
#' @param options Options for the layout generator, a named list. See details
#' below.
#' @param weights Optional edge weights. Supply \code{NULL} here if you want to
#' weight edges equally. By default the \code{weight} edge attribute is used if
#' the graph has one.
#' @param fixed Logical vector, it can be used to fix some vertices. All
#' vertices for which it is \code{TRUE} are kept at the coordinates supplied in
#' the \code{seed} matrix. It is ignored it \code{NULL} or if \code{use.seed}
#' is \code{FALSE}.
#' @param dim Either \sQuote{2} or \sQuote{3}, it specifies whether we want a
#' two dimensional or a three dimensional layout. Note that because of the
#' nature of the DrL algorithm, the three dimensional layout takes
#' significantly longer to compute.
#' @return A numeric matrix with two columns.
#' @author Shawn Martin (\url{http://www.cs.otago.ac.nz/homepages/smartin/})
#' and Gabor Csardi \email{csardi.gabor@@gmail.com} for the R/igraph interface
#' and the three dimensional version.
#' @seealso \code{\link{layout}} for other layout generators.
#' @references See the following technical report: Martin, S., Brown, W.M.,
#' Klavans, R., Boyack, K.W., DrL: Distributed Recursive (Graph) Layout. SAND
#' Reports, 2008. 2936: p. 1-10.
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- as.undirected(sample_pa(100, m=1))
#' l <- layout_with_drl(g, options=list(simmer.attraction=0))
#' \dontrun{
#' plot(g, layout=l, vertex.size=3, vertex.label=NA)
#' }
#' 
layout_with_drl <- function(graph, use.seed = FALSE,
                       seed=matrix(runif(vcount(graph)*2), ncol=2),
                       options=drl_defaults$default,
                       weights=E(graph)$weight,
                       fixed=NULL,
                       dim=2)
{
    if (!is_igraph(graph)) {
        stop("Not a graph object")
    }
    if (dim != 2 && dim != 3) {
      stop("`dim' must be 2 or 3")
    }
    use.seed <- as.logical(use.seed)
    seed <- as.matrix(seed)
    options.tmp <- drl_defaults$default
    options.tmp[names(options)] <- options
    options <- options.tmp
    if (!is.null(weights)) {
      weights <- as.numeric(weights)
    }
    if (!is.null(fixed)) {
      fixed <- as.logical(fixed)
    }
    on.exit(.Call("R_igraph_finalizer", PACKAGE = "igraph"))
    if (dim==2) {
      res <- .Call("R_igraph_layout_drl", graph, seed, use.seed, options, 
                   weights, fixed, PACKAGE = "igraph")
    } else {
      res <- .Call("R_igraph_layout_drl_3d", graph, seed, use.seed, options, 
                   weights, fixed, PACKAGE = "igraph")
    }      
    res
}


#' @rdname layout_with_drl
#' @param ... Passed to \code{layout_with_drl}.
#' @export

with_drl <- function(...) layout_spec(layout_with_drl, ...)


#' @export

igraph.drl.default <- list(edge.cut=32/40,
                           init.iterations=0,
                           init.temperature=2000,
                           init.attraction=10,
                           init.damping.mult=1.0,
                           liquid.iterations=200,
                           liquid.temperature=2000,
                           liquid.attraction=10,
                           liquid.damping.mult=1.0,
                           expansion.iterations=200,
                           expansion.temperature=2000,
                           expansion.attraction=2,
                           expansion.damping.mult=1.0,
                           cooldown.iterations=200,
                           cooldown.temperature=2000,
                           cooldown.attraction=1,
                           cooldown.damping.mult=.1,
                           crunch.iterations=50,
                           crunch.temperature=250,
                           crunch.attraction=1,
                           crunch.damping.mult=0.25,
                           simmer.iterations=100,
                           simmer.temperature=250,
                           simmer.attraction=.5,
                           simmer.damping.mult=0)

#' @export

igraph.drl.coarsen <- list(edge.cut=32/40,
                           init.iterations=0,
                           init.temperature=2000,
                           init.attraction=10,
                           init.damping.mult=1.0,
                           liquid.iterations=200,
                           liquid.temperature=2000,
                           liquid.attraction=2,
                           liquid.damping.mult=1.0,
                           expansion.iterations=200,
                           expansion.temperature=2000,
                           expansion.attraction=10,
                           expansion.damping.mult=1.0,
                           cooldown.iterations=200,
                           cooldown.temperature=2000,
                           cooldown.attraction=1,
                           cooldown.damping.mult=.1,
                           crunch.iterations=50,
                           crunch.temperature=250,
                           crunch.attraction=1,
                           crunch.damping.mult=0.25,
                           simmer.iterations=100,
                           simmer.temperature=250,
                           simmer.attraction=.5,
                           simmer.damping.mult=0)

#' @export

igraph.drl.coarsest <- list(edge.cut=32/40,
                            init.iterations=0,
                            init.temperature=2000,
                            init.attraction=10,
                            init.damping.mult=1.0,
                            liquid.iterations=200,
                            liquid.temperature=2000,
                            liquid.attraction=2,
                            liquid.damping.mult=1.0,
                            expansion.iterations=200,
                            expansion.temperature=2000,
                            expansion.attraction=10,
                            expansion.damping.mult=1.0,
                            cooldown.iterations=200,
                            cooldown.temperature=2000,
                            cooldown.attraction=1,
                            cooldown.damping.mult=.1,
                            crunch.iterations=200,
                            crunch.temperature=250,
                            crunch.attraction=1,
                            crunch.damping.mult=0.25,
                            simmer.iterations=100,
                            simmer.temperature=250,
                            simmer.attraction=.5,
                            simmer.damping.mult=0)

#' @export

igraph.drl.refine <- list(edge.cut=32/40,
                          init.iterations=0,
                          init.temperature=50,
                          init.attraction=.5,
                          init.damping.mult=1.0,
                          liquid.iterations=0,
                          liquid.temperature=2000,
                          liquid.attraction=2,
                          liquid.damping.mult=1.0,
                          expansion.iterations=50,
                          expansion.temperature=500,
                          expansion.attraction=.1,
                          expansion.damping.mult=.25,
                          cooldown.iterations=50,
                          cooldown.temperature=250,
                          cooldown.attraction=1,
                          cooldown.damping.mult=.1,
                          crunch.iterations=50,
                          crunch.temperature=250,
                          crunch.attraction=1,
                          crunch.damping.mult=0.25,
                          simmer.iterations=0,
                          simmer.temperature=250,
                          simmer.attraction=.5,
                          simmer.damping.mult=0)

#' @export

igraph.drl.final <- list(edge.cut=32/40,
                         init.iterations=0,
                         init.temperature=50,
                         init.attraction=.5,
                         init.damping.mult=0,
                         liquid.iterations=0,
                         liquid.temperature=2000,
                         liquid.attraction=2,
                         liquid.damping.mult=1.0,
                         expansion.iterations=50,
                         expansion.temperature=2000,
                         expansion.attraction=2,
                         expansion.damping.mult=1.0,
                         cooldown.iterations=50,
                         cooldown.temperature=200,
                         cooldown.attraction=1,
                         cooldown.damping.mult=.1,
                         crunch.iterations=50,
                         crunch.temperature=250,
                         crunch.attraction=1,
                         crunch.damping.mult=0.25,
                         simmer.iterations=25,
                         simmer.temperature=250,
                         simmer.attraction=.5,
                         simmer.damping.mult=0)

#' @export

drl_defaults <- list(
  coarsen = igraph.drl.coarsen,
  coarsest = igraph.drl.coarsest,
  default = igraph.drl.default,
  final = igraph.drl.final,
  refine = igraph.drl.refine
)
