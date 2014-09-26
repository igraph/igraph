#   IGraph R package
#   Copyright (C) 2003-2012  Gabor Csardi <csardi.gabor@gmail.com>
#   334 Harvard street, Cambridge, MA 02139 USA
#   
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#   
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
#   02110-1301 USA
#
###################################################################

###################################################################
# Layouts
###################################################################



#' Randomly place vertices on a plane or in 3d space
#' 
#' This function uniformly randomly places the vertices of the graph in two or
#' three dimensions.
#' 
#' Randomly places vertices on a [-1,1] square (in 2d) or in a cube (in 3d). It
#' is probably a useless layout, but it can use as a starting point for other
#' layout generators.
#' 
#' @param graph The input graph.
#' @param dim Integer scalar, the dimension of the space to use. It must be 2
#' or 3.
#' @return A numeric matrix with two or three columns.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @keywords graphs
layout.random <- function(graph, dim=2) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (dim==2) {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    .Call("R_igraph_layout_random", graph,
          PACKAGE="igraph")
  } else if (dim==3) {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    .Call("R_igraph_layout_random_3d", graph,
          PACKAGE="igraph")
  } else {
    stop("Invalid `dim' value");
  }
}



#' Graph layout with vertices on a circle.
#' 
#' Place vertices on a circle, in the order of their vertex ids.
#' 
#' If you want to order the vertices differently, then permute them using the
#' \code{\link{permute.vertices}} function.
#' 
#' @param graph The input graph.
#' @param order The vertices to place on the circle, in the order of their
#' desired placement. Vertices that are not included here will be placed at
#' (0,0).
#' @return A numeric matrix with two columns, and one row for each vertex.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @keywords graphs
#' @examples
#' 
#' ## Place vertices on a circle, order them according to their
#' ## community
#' library(igraphdata)
#' data(karate)
#' karate_groups <- optimal.community(karate)
#' coords <- layout.circle(karate, order =
#'           order(membership(karate_groups)))
#' V(karate)$label <- sub("Actor ", "", V(karate)$name)
#' V(karate)$label.color <- membership(karate_groups)
#' V(karate)$shape <- "none"
#' plot(karate, layout = coords)
#' 
layout.circle <- function(graph, order=V(graph)) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  order <- as.igraph.vs(graph, order) - 1L
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_layout_circle", graph, order,
        PACKAGE="igraph")
}



#' Graph layout with vertices on the surface of a sphere
#' 
#' Place vertices on a sphere, approximately uniformly, in the order of their
#' vertex ids.
#' 
#' \code{layout.sphere} places the vertices (approximately) uniformly on the
#' surface of a sphere, this is thus a 3d layout. It is not clear however what
#' \dQuote{uniformly on a sphere} means.
#' 
#' If you want to order the vertices differently, then permute them using the
#' \code{\link{permute.vertices}} function.
#' 
#' @param graph The input graph.
#' @return A numeric matrix with three columns, and one row for each vertex.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @keywords graphs
layout.sphere <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_layout_sphere", graph,
        PACKAGE="igraph")
}  



#' The graphopt layout algorithm
#' 
#' A force-directed layout algorithm, that scales relatively well to large
#' graphs.
#' 
#' \code{layout.graphopt} is a port of the graphopt layout algorithm by Michael
#' Schmuhl. graphopt version 0.4.1 was rewritten in C and the support for
#' layers was removed (might be added later) and a code was a bit reorganized
#' to avoid some unneccessary steps is the node charge (see below) is zero.
#' 
#' graphopt uses physical analogies for defining attracting and repelling
#' forces among the vertices and then the physical system is simulated until it
#' reaches an equilibrium. (There is no simulated annealing or anything like
#' that, so a stable fixed point is not guaranteed.)
#' 
#' See also \url{http://www.schmuhl.org/graphopt/} for the original graphopt.
#' 
#' @param graph The input graph.
#' @param start If given, then it should be a matrix with two columns and one
#' line for each vertex. This matrix will be used as starting positions for the
#' algorithm. If not given, then a random starting matrix is used.
#' @param niter Integer scalar, the number of iterations to perform.  Should be
#' a couple of hundred in general. If you have a large graph then you might
#' want to only do a few iterations and then check the result. If it is not
#' good enough you can feed it in again in the \code{start} argument. The
#' default value is 500.
#' @param charge The charge of the vertices, used to calculate electric
#' repulsion. The default is 0.001.
#' @param mass The mass of the vertices, used for the spring forces. The
#' default is 30.
#' @param spring.length The length of the springs, an integer number. The
#' default value is zero.
#' @param spring.constant The spring constant, the default value is one.
#' @param max.sa.movement Real constant, it gives the maximum amount of
#' movement allowed in a single step along a single axis. The default value is
#' 5.
#' @return A numeric matrix with two columns, and a row for each vertex.
#' @author Michael Schmuhl for the original graphopt code, rewritten and
#' wrapped by Gabor Csardi \email{csardi.gabor@@gmail.com}.
#' @keywords graphs
layout.graphopt <- function(graph, start=NULL, niter=500, charge=0.001,
                            mass=30, spring.length=0, spring.constant=1,
                            max.sa.movement=5) {
  
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (!is.null(start)) {
    start <- structure(as.numeric(start), dim=dim(start))
  }
  niter <- as.double(niter)
  charge <- as.double(charge)
  mass <- as.double(mass)
  spring.length <- as.double(spring.length)
  spring.constant <- as.double(spring.constant)
  max.sa.movement <- as.double(max.sa.movement)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_layout_graphopt", graph, niter, charge, mass,
        spring.length, spring.constant, max.sa.movement, start,
        PACKAGE="igraph")
}



#' Large Graph Layout
#' 
#' A layout generator for larger graphs.
#' 
#' \code{layout.lgl} is for large connected graphs, it is similar to the layout
#' generator of the Large Graph Layout software
#' (\url{http://lgl.sourceforge.net/}).
#' 
#' @param graph The input graph
#' @param maxiter The maximum number of iterations to perform (150).
#' @param maxdelta The maximum change for a vertex during an iteration (the
#' number of vertices).
#' @param area The area of the surface on which the vertices are placed (square
#' of the number of vertices).
#' @param coolexp The cooling exponent of the simulated annealing (1.5).
#' @param repulserad Cancellation radius for the repulsion (the \code{area}
#' times the number of vertices).
#' @param cellsize The size of the cells for the grid. When calculating the
#' repulsion forces between vertices only vertices in the same or neighboring
#' grid cells are taken into account (the fourth root of the number of
#' \code{area}.
#' @param root The id of the vertex to place at the middle of the layout. The
#' default value is -1 which means that a random vertex is selected.
#' @return A numeric matrix with two columns and as many rows as vertices.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @keywords graphs
layout.lgl <- function(graph, maxiter=150, maxdelta=vcount(graph),
                       area=vcount(graph)^2, coolexp=1.5,
                       repulserad=area * vcount(graph),
                       cellsize=sqrt(sqrt(area)), root=NULL) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (is.null(root)) {
    root <- -1
  } else {
    root <- as.igraph.vs(graph, root)-1
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_layout_lgl", graph, as.double(maxiter),
        as.double(maxdelta), as.double(area), as.double(coolexp),
        as.double(repulserad), as.double(cellsize), root,
        PACKAGE="igraph")
}



#' The Reingold-Tilford graph layout algorithm
#' 
#' A tree-like layout, it is perfect for trees, acceptable for graphs with not
#' too many cycles.
#' 
#' Arranges the nodes in a tree where the given node is used as the root.  The
#' tree is directed downwards and the parents are centered above its children.
#' For the exact algorithm, the refernce below.
#' 
#' If the given graph is not a tree, a breadth-first search is executed first
#' to obtain a possible spanning tree.
#' 
#' @param graph The input graph.
#' @param root The index of the root vertex or root vertices.  If this is a
#' non-empty vector then the supplied vertex ids are used as the roots of the
#' trees (or a single tree if the graph is connected).  If it is an empty
#' vector, then the root vertices are automatically calculated based on
#' topological sorting, performed with the opposite mode than the \code{mode}
#' argument. After the vertices have been sorted, one is selected from each
#' component.
#' @param circular Logical scalar, whether to plot the tree in a circular
#' fashion. Defaults to \code{FALSE}, so the tree branches are going bottom-up
#' (or top-down, see the \code{flip.y} argument.
#' @param rootlevel This argument can be useful when drawing forests which are
#' not trees (i.e. they are unconnected and have tree components). It specifies
#' the level of the root vertices for every tree in the forest. It is only
#' considered if the \code{roots} argument is not an empty vector.
#' @param mode Specifies which edges to consider when building the tree.  If it
#' is \sQuote{out}, then only the outgoing, if it is \sQuote{in}, then only the
#' incoming edges of a parent are considered. If it is \sQuote{all} then all
#' edges are used (this was the behavior in igraph 0.5 and before). This
#' parameter also influences how the root vertices are calculated, if they are
#' not given. See the \code{roots} parameter.
#' @param flip.y Logical scalar, whether to flip the \sQuote{y} coordinates.
#' The default is flipping because that puts the root vertex on the top.
#' @return A numeric matrix with two columns, and one row for each vertex.
#' @author Tamas Nepusz \email{ntamas@@gmail.com} and Gabor Csardi
#' \email{csardi.gabor@@gmail.com}
#' @references Reingold, E and Tilford, J (1981). Tidier drawing of trees.
#' \emph{IEEE Trans. on Softw. Eng.}, SE-7(2):223--228.
#' @keywords graphs
#' @examples
#' 
#' tree <- graph.tree(20, 3)
#' plot(tree, layout=layout.reingold.tilford)
#' plot(tree, layout=layout.reingold.tilford(tree, flip.y=FALSE))
#' plot(tree, layout=layout.reingold.tilford(tree, circular=TRUE))
#' 
#' tree2 <- graph.tree(10, 3) + graph.tree(10, 2)
#' plot(tree2, layout=layout.reingold.tilford)
#' plot(tree2, layout=layout.reingold.tilford(tree2, root=c(1,11),
#'                                            rootlevel=c(2,1)))
#' 
layout.reingold.tilford <- function(graph, root=numeric(), circular=FALSE,
                                    rootlevel=numeric(), mode="out",
                                    flip.y=TRUE) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  root <- as.igraph.vs(graph, root)-1
  circular <- as.logical(circular)
  rootlevel <- as.double(rootlevel)
  mode <- switch(igraph.match.arg(mode), "out"=1, "in"=2, "all"=3,
                 "total"=3)
  flip.y <- as.logical(flip.y)
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_layout_reingold_tilford", graph, root, mode,
               rootlevel, circular, PACKAGE="igraph")
  if (flip.y) { res[,2] <- max(res[,2])-res[,2] }
  res
}



#' Merging graph layouts
#' 
#' Place several graphs on the same layout
#' 
#' \code{layout.merge} takes a list of graphs and a list of coordinates and
#' places the graphs in a common layout. The method to use is chosen via the
#' \code{method} parameter, although right now only the \code{dla} method is
#' implemented.
#' 
#' The \code{dla} method covers the graph with circles.  Then it sorts the
#' graphs based on the number of vertices first and places the largest graph at
#' the center of the layout. Then the other graphs are placed in decreasing
#' order via a DLA (diffision limited aggregation) algorithm: the graph is
#' placed randomly on a circle far away from the center and a random walk is
#' conducted until the graph walks into the larger graphs already placed or
#' walks too far from the center of the layout.
#' 
#' The \code{piecewise.layout} function disassembles the graph first into
#' maximal connected components and calls the supplied \code{layout} function
#' for each component separately. Finally it merges the layouts via calling
#' \code{layout.merge}.
#' 
#' @aliases layout.merge piecewise.layout
#' @param graphs A list of graph objects.
#' @param layouts A list of two-column matrices.
#' @param method Character constant giving the method to use. Right now only
#' \code{dla} is implemented.
#' @param graph The input graph.
#' @param layout A function object, the layout function to use.
#' @param \dots Additional arguments to pass to the \code{layout} layout
#' function.
#' @return A matrix with two columns and as many lines as the total number of
#' vertices in the graphs.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{plot.igraph}}, \code{\link{tkplot}},
#' \code{\link{layout}}, \code{\link{graph.disjoint.union}}
#' @keywords graphs
#' @examples
#' 
#' # create 20 scale-free graphs and place them in a common layout
#' graphs <- lapply(sample(5:20, 20, replace=TRUE),
#'           barabasi.game, directed=FALSE)
#' layouts <- lapply(graphs, layout.kamada.kawai)
#' lay <- layout.merge(graphs, layouts)
#' g <- graph.disjoint.union(graphs)
#' \dontrun{plot(g, layout=lay, vertex.size=3, labels=NA, edge.color="black")}
#' 
layout.merge <- function(graphs, layouts, method="dla") {

  if (!all(sapply(graphs, is.igraph))) {
    stop("Not a graph object")
  }
  if (method == "dla") {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    res <- .Call("R_igraph_layout_merge_dla",
                 graphs, layouts,
                 PACKAGE="igraph")
  } else {
    stop("Invalid `method'.")
  }
  res
}



#' Normalize coordinates for plotting graphs
#' 
#' Rescale coordinates linearly to be within given bounds.
#' 
#' \code{layout.norm} normalizes a layout, it linearly transforms each
#' coordinate separately to fit into the given limits.
#' 
#' @param layout A matrix with two or three columns, the layout to normalize.
#' @param xmin,xmax The limits for the first coordinate, if one of them or both
#' are \code{NULL} then no normalization is performed along this direction.
#' @param ymin,ymax The limits for the second coordinate, if one of them or
#' both are \code{NULL} then no normalization is performed along this
#' direction.
#' @param zmin,zmax The limits for the third coordinate, if one of them or both
#' are \code{NULL} then no normalization is performed along this direction.
#' @return A numeric matrix with at the same dimension as \code{layout}.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @keywords graphs
layout.norm <- function(layout, xmin=-1, xmax=1, ymin=-1, ymax=1,
                          zmin=-1, zmax=1) {

  if (!is.matrix(layout)) {
    stop("`layout' not a matrix")
  }
  if (ncol(layout) != 2 && ncol(layout) != 3) {
    stop("`layout' should have 2 or three columns")
  }
  
  if (!is.null(xmin) && !is.null(xmax)) {
    layout[,1] <- .layout.norm.col(layout[,1], xmin, xmax)
  }

  if (!is.null(ymin) && !is.null(ymax)) {
    layout[,2] <- .layout.norm.col(layout[,2], ymin, ymax)
  }
  
  if (ncol(layout)==3 && !is.null(zmin) && !is.null(zmax)) {
    layout[,3] <- .layout.norm.col(layout[,3], zmin, zmax)
  }

  layout
}

.layout.norm.col <- function(v, min, max) {

  vr <- range(v)
  if (vr[1]==vr[2]) {
    fac <- 1
  } else {
    fac <- (max-min)/(vr[2]-vr[1])
  }

  (v-vr[1]) * fac + min
}

piecewise.layout <- function(graph, layout=layout.kamada.kawai, ...) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  
  V(graph)$id <- seq(vcount(graph))
  gl <- decompose.graph(graph)
  ll <- lapply(gl, layout, ...)
  
  l <- layout.merge(gl, ll)
  l[ unlist(sapply(gl, get.vertex.attribute, "id")), ] <- l[]
  l
}



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
#' \code{igraph.drl.default}, \code{igraph.drl.coarsen},
#' \code{igraph.drl.coarsest}, \code{igraph.drl.refine} and
#' \code{igraph.drl.final}.  }
#' 
#' @aliases layout.drl igraph.drl.default igraph.drl.coarsen
#' igraph.drl.coarsest igraph.drl.refine igraph.drl.final
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
#' @keywords graphs
#' @examples
#' 
#' g <- as.undirected(ba.game(100, m=1))
#' l <- layout.drl(g, options=list(simmer.attraction=0))
#' \dontrun{
#' plot(g, layout=l, vertex.size=3, vertex.label=NA)
#' }
#' 
layout.drl <- function(graph, use.seed = FALSE,
                       seed=matrix(runif(vcount(graph)*2), ncol=2),
                       options=igraph.drl.default,
                       weights=E(graph)$weight,
                       fixed=NULL,
                       dim=2)
{
    if (!is.igraph(graph)) {
        stop("Not a graph object")
    }
    if (dim != 2 && dim != 3) {
      stop("`dim' must be 2 or 3")
    }
    use.seed <- as.logical(use.seed)
    seed <- as.matrix(seed)
    options.tmp <- igraph.drl.default
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



#' Choose an appropriate graph layout algorithm automaticall
#' 
#' This function tries to choose an appropriate graph layout algorithm for the
#' graph, automatically, based on a simple algorithm. See details below.
#' 
#' \code{layout.auto} tries to choose an appropriate layout function for the
#' supplied graph, and uses that to generate the layout. The current
#' implementation works like this: \enumerate{ \item If the graph has a graph
#' attribute called \sQuote{layout}, then this is used. If this attribute is an
#' R function, then it is called, with the graph and any other extra arguments.
#' \item Otherwise, if the graph has vertex attributes called \sQuote{x} and
#' \sQuote{y}, then these are used as coordinates. If the graph has an
#' additional \sQuote{z} vertex attribute, that is also used.  \item Otherwise,
#' if the graph is connected and has less than 1000 vertices, the Kamada-Kawai
#' layout is used, by calling \code{layout.kamada.kawai}.  \item Otherwise the
#' DrL layout is used, \code{layout.drl} is called.  }
#' 
#' @param graph The input graph
#' @param dim Dimensions, should be 2 or 3.
#' @param \dots The extra arguments are passed to the real layout function.
#' @return A numeric matrix with two or three columns.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{plot.igraph}}
#' @keywords graphs
layout.auto <- function(graph, dim=2, ...) {

  ## 1. If there is a 'layout' graph attribute, we just use that.
  ## 2. Otherwise, if there are vertex attributes called 'x' and 'y',
  ##    we use those (and the 'z' vertex attribute as well, if present).
  ## 3. Otherwise, if the graph is small (<1000) we use
  ##    the Kamada-Kawai layout.
  ## 5. Otherwise we use the DrL layout generator.
  
  if ("layout" %in% list.graph.attributes(graph)) {
    lay <- get.graph.attribute(graph, "layout")
    if (is.function(lay)) {
      lay(graph, ...)
    } else {
      lay
    }

  } else if ( all(c("x", "y") %in% list.vertex.attributes(graph)) ) {
    if ("z" %in% list.vertex.attributes(graph)) {
      cbind(V(graph)$x, V(graph)$y, V(graph)$z)
    } else {
      cbind(V(graph)$x, V(graph)$y)
    }

  } else if (vcount(graph) < 1000) {
    layout.kamada.kawai(graph, dim=dim, ...)

  } else {
    layout.drl(graph, dim=dim, ...)
  }
  
}



#' The Sugiyama graph layout generator
#' 
#' Sugiyama layout algorithm for layered directed acyclic graphs. The algorithm
#' minimized edge crossings.
#' 
#' This layout algorithm is designed for directed acyclic graphs where each
#' vertex is assigned to a layer. Layers are indexed from zero, and vertices of
#' the same layer will be placed on the same horizontal line. The X coordinates
#' of vertices within each layer are decided by the heuristic proposed by
#' Sugiyama et al. to minimize edge crossings.
#' 
#' You can also try to lay out undirected graphs, graphs containing cycles, or
#' graphs without an a priori layered assignment with this algorithm. igraph
#' will try to eliminate cycles and assign vertices to layers, but there is no
#' guarantee on the quality of the layout in such cases.
#' 
#' The Sugiyama layout may introduce \dQuote{bends} on the edges in order to
#' obtain a visually more pleasing layout. This is achieved by adding dummy
#' nodes to edges spanning more than one layer. The resulting layout assigns
#' coordinates not only to the nodes of the original graph but also to the
#' dummy nodes. The layout algorithm will also return the extended graph with
#' the dummy nodes.
#' 
#' For more details, see the reference below.
#' 
#' @param graph The input graph.
#' @param layers A numeric vector or \code{NULL}. If not \code{NULL}, then it
#' should specify the layer index of the vertices. Layers are numbered from
#' one. If \code{NULL}, then igraph calculates the layers automatically.
#' @param hgap Real scalar, the minimum horizontal gap between vertices in the
#' same layer.
#' @param vgap Real scalar, the distance between layers.
#' @param maxiter Integer scalar, the maximum number of iterations in the
#' crossing minimization stage. 100 is a reasonable default; if you feel that
#' you have too many edge crossings, increase this.
#' @param weights Optional edge weight vector. If \code{NULL}, then the
#' 'weight' edge attribute is used, if there is one. Supply \code{NA} here and
#' igraph ignores the edge weights.
#' @param attributes Which graph/vertex/edge attributes to keep in the extended
#' graph. \sQuote{default} keeps the \sQuote{size}, \sQuote{size2},
#' \sQuote{shape}, \sQuote{label} and \sQuote{color} vertex attributes and the
#' \sQuote{arrow.mode} and \sQuote{arrow.size} edge attributes. \sQuote{all}
#' keep all graph, vertex and edge attributes, \sQuote{none} keeps none of
#' them.
#' @return A list with the components: \item{layout}{The layout, a two-column
#' matrix, for the original graph vertices.} \item{layout.dummy}{The layout for
#' the dummy vertices, a two column matrix.} \item{extd_graph}{The original
#' graph, extended with dummy vertices.  The \sQuote{dummy} vertex attribute is
#' set on this graph, it is a logical attributes, and it tells you whether the
#' vertex is a dummy vertex. The \sQuote{layout} graph attribute is also set,
#' and it is the layout matrix for all (original and dummy) vertices.}
#' @author Tamas Nepusz \email{ntamas@@gmail.com}
#' @references K. Sugiyama, S. Tagawa and M. Toda, "Methods for Visual
#' Understanding of Hierarchical Systems". IEEE Transactions on Systems, Man
#' and Cybernetics 11(2):109-125, 1981.
#' @keywords graphs
#' @examples
#' 
#' ## Data taken from http://tehnick-8.narod.ru/dc_clients/
#' DC <- graph.formula("DC++" -+
#'                 "LinuxDC++":"BCDC++":"EiskaltDC++":"StrongDC++":"DiCe!++",
#'                 "LinuxDC++" -+ "FreeDC++", "BCDC++" -+ "StrongDC++",
#'                 "FreeDC++" -+ "BMDC++":"EiskaltDC++",
#'                 "StrongDC++" -+ "AirDC++":"zK++":"ApexDC++":"TkDC++",
#'                 "StrongDC++" -+ "StrongDC++ SQLite":"RSX++",
#'                 "ApexDC++" -+ "FlylinkDC++ ver <= 4xx",
#'                 "ApexDC++" -+ "ApexDC++ Speed-Mod":"DiCe!++",
#'                 "StrongDC++ SQLite" -+ "FlylinkDC++ ver >= 5xx",
#'                 "ApexDC++ Speed-Mod" -+ "FlylinkDC++ ver <= 4xx",
#'                 "ApexDC++ Speed-Mod" -+ "GreylinkDC++",
#'                 "FlylinkDC++ ver <= 4xx" -+ "FlylinkDC++ ver >= 5xx",
#'                 "FlylinkDC++ ver <= 4xx" -+ AvaLink,
#'                 "GreylinkDC++" -+ AvaLink:"RayLinkDC++":"SparkDC++":PeLink)
#' 
#' ## Use edge types
#' E(DC)$lty <- 1
#' E(DC)["BCDC++" %->% "StrongDC++"]$lty <- 2
#' E(DC)["FreeDC++" %->% "EiskaltDC++"]$lty <- 2
#' E(DC)["ApexDC++" %->% "FlylinkDC++ ver <= 4xx"]$lty <- 2
#' E(DC)["ApexDC++" %->% "DiCe!++"]$lty <- 2
#' E(DC)["StrongDC++ SQLite" %->% "FlylinkDC++ ver >= 5xx"]$lty <- 2
#' E(DC)["GreylinkDC++" %->% "AvaLink"]$lty <- 2
#' 
#' ## Layers, as on the plot
#' layers <- list(c("DC++"),
#'                c("LinuxDC++", "BCDC++"),
#'                c("FreeDC++", "StrongDC++"),
#'                c("BMDC++", "EiskaltDC++", "AirDC++", "zK++", "ApexDC++",
#'                  "TkDC++", "RSX++"),
#'                c("StrongDC++ SQLite", "ApexDC++ Speed-Mod", "DiCe!++"),
#'                c("FlylinkDC++ ver <= 4xx", "GreylinkDC++"),
#'                c("FlylinkDC++ ver >= 5xx", "AvaLink", "RayLinkDC++",
#'                  "SparkDC++", "PeLink"))
#' 
#' ## Check that we have all nodes
#' all(sort(unlist(layers)) == sort(V(DC)$name))
#' 
#' ## Add some graphical parameters
#' V(DC)$color <- "white"
#' V(DC)$shape <- "rectangle"
#' V(DC)$size <- 20
#' V(DC)$size2 <- 10
#' V(DC)$label <- lapply(V(DC)$name, function(x)
#'                       paste(strwrap(x, 12), collapse="\n"))
#' E(DC)$arrow.size <- 0.5
#' 
#' ## Create a similar layout using the predefined layers
#' lay1 <- layout.sugiyama(DC, layers=apply(sapply(layers,
#'                         function(x) V(DC)$name %in% x), 1, which))
#' 
#' ## Simple plot, not very nice
#' par(mar=rep(.1, 4))
#' plot(DC, layout=lay1$layout, vertex.label.cex=0.5)
#' 
#' ## Sugiyama plot
#' plot(lay1$extd_graph, vertex.label.cex=0.5)
#' 
#' ## The same with automatic layer calculation
#' ## Keep vertex/edge attributes in the extended graph
#' lay2 <- layout.sugiyama(DC, attributes="all")
#' plot(lay2$extd_graph, vertex.label.cex=0.5)
#' 
#' ## Another example, from the following paper:
#' ## Markus Eiglsperger, Martin Siebenhaller, Michael Kaufmann:
#' ## An Efficient Implementation of Sugiyama's Algorithm for
#' ## Layered Graph Drawing, Journal of Graph Algorithms and
#' ## Applications 9, 305--325 (2005).
#' 
#' ex <- graph.formula( 0 -+ 29: 6: 5:20: 4,
#'                      1 -+ 12,
#'                      2 -+ 23: 8,
#'                      3 -+  4,
#'                      4,
#'                      5 -+  2:10:14:26: 4: 3,
#'                      6 -+  9:29:25:21:13,
#'                      7,
#'                      8 -+ 20:16,
#'                      9 -+ 28: 4,
#'                     10 -+ 27,
#'                     11 -+  9:16,
#'                     12 -+  9:19,
#'                     13 -+ 20,
#'                     14 -+ 10,
#'                     15 -+ 16:27,
#'                     16 -+ 27,
#'                     17 -+  3,
#'                     18 -+ 13,
#'                     19 -+  9,
#'                     20 -+  4,
#'                     21 -+ 22,
#'                     22 -+  8: 9,
#'                     23 -+  9:24,
#'                     24 -+ 12:15:28,
#'                     25 -+ 11,
#'                     26 -+ 18,
#'                     27 -+ 13:19,
#'                     28 -+  7,
#'                     29 -+ 25                    )
#' 
#' layers <- list( 0, c(5, 17), c(2, 14, 26, 3), c(23, 10, 18), c(1, 24),
#'                 12, 6, c(29,21), c(25,22), c(11,8,15), 16, 27, c(13,19),
#'                 c(9, 20), c(4, 28), 7 )
#' 
#' layex <- layout.sugiyama(ex, layers=apply(sapply(layers,
#'                         function(x) V(ex)$name %in% as.character(x)),
#'                         1, which))
#' 
#' origvert <- c(rep(TRUE, vcount(ex)), rep(FALSE, nrow(layex$layout.dummy)))
#' realedge <- get.edgelist(layex$extd_graph)[,2] <= vcount(ex)
#' plot(layex$extd_graph, vertex.label.cex=0.5,
#'      edge.arrow.size=.5, 
#'      vertex.size=ifelse(origvert, 5, 0),
#'      vertex.shape=ifelse(origvert, "square", "none"),
#'      vertex.label=ifelse(origvert, V(ex)$name, ""),
#'      edge.arrow.mode=ifelse(realedge, 2, 0))
#' 
layout.sugiyama <- function(graph, layers=NULL, hgap=1, vgap=1,
                            maxiter=100, weights=NULL,


#' Graph, vertex and edge attributes
#' 
#' Attributes are associated values belonging to a graph, vertices or edges.
#' These can represent some property, like data about how the graph was
#' constructed, the color of the vertices when the graph is plotted, or simply
#' the weights of the edges in a weighted graph.
#' 
#' There are three types of attributes in igraph: graph, vertex and edge
#' attributes. Graph attributes are associated with graph, vertex attributes
#' with vertices and edge attributes with edges.
#' 
#' Examples for graph attributes are the date when the graph data was collected
#' or other types of memos like the type of the data, or whether the graph is a
#' simple graph, ie. one without loops and multiple edges.
#' 
#' Examples of vertex attributes are vertex properties, like the vertex
#' coordinates for the visualization of the graph, or other visualization
#' parameters, or meta-data associated with the vertices, like the gender and
#' the age of the individuals in a friendship network, the type of the neurons
#' in a graph representing neural circuitry or even some pre-computed structual
#' properties, like the betweenness centrality of the vertices.
#' 
#' Examples of edge attributes are data associated with edges: most commonly
#' edge weights, or visualization parameters.
#' 
#' In recent igraph versions, arbitrary R objects can be assigned as graph,
#' vertex or edge attributes.
#' 
#' Some igraph functions use the values or graph, vertex and edge attributes if
#' they are present but this is not done in the current version very
#' extensively. Expect more in the (near) future.
#' 
#' Graph attributes can be created with the \code{set.graph.attribute}
#' function, and removed with \code{remove.graph.attribute}. Graph attributes
#' are queried with \code{get.graph.attribute} and the assigned graph
#' attributes are listed with \code{list.graph.attributes}.
#' 
#' The function \code{graph.attributes} returns all graph attributes in a named
#' list. This function has a counterpart that sets all graph attributes at
#' once, see an example below.
#' 
#' There is a simpler notation for using graph attributes: the
#' \sQuote{\code{$}} operator. It can be used both to query and set graph
#' attributes, see an example below.
#' 
#' The functions for vertex attributes are \code{set.vertex.attribute},
#' \code{get.vertex.attribute}, \code{remove.vertex.attribute} and
#' \code{list.vertex.attributes} and for edge attributes they are
#' \code{set.edge.attribute}, \code{get.edge.attribute},
#' \code{remove.edge.attribute} and \code{list.edge.attributes}.
#' 
#' Similarly to graph attributes, \code{vertex.attributes} returns all vertex
#' attributes, in a named list, and \code{edge.attributes} returns all edge
#' attributes, in a named list.
#' 
#' There is however a (syntactically) much simpler way to handle vertex and
#' edge attribute by using vertex and edge selectors, it works like this:
#' \code{V(g)} selects all vertices in a graph, and \code{V(g)$name} queries
#' the \code{name} attribute for all vertices. Similarly is \code{vs} is a
#' vertex set \code{vs$name} gives the values of the \code{name} attribute for
#' the vertices in the vertex set.
#' 
#' This form can also be used to set the values of the attributes, like the
#' regular R convention: \preformatted{V(g)$color <- "red"} It works for vertex
#' subsets as well: \preformatted{V(g)[1:5]$color <- "green"}
#' 
#' The notation for edges is similar: \code{E(g)} means all edges
#' \code{E(g)$weight} is the \code{weight} attribute for all edges, etc.
#' 
#' See also the manual page for \code{iterators} about how to create various
#' vertex and edge sets.
#' 
#' @aliases attributes set.graph.attribute get.graph.attribute
#' remove.graph.attribute list.graph.attributes graph.attributes
#' graph.attributes<- set.vertex.attribute get.vertex.attribute
#' remove.vertex.attribute list.vertex.attributes vertex.attributes
#' vertex.attributes<- set.edge.attribute get.edge.attribute
#' remove.edge.attribute list.edge.attributes edge.attributes edge.attributes<-
#' @param graph The graph object to work on. Note that the original graph is
#' never modified, a new graph object is returned instead; if you don't assign
#' it to a variable your modifications will be lost! See examples below.
#' @param name Character constant, the name of the attribute.
#' @param index Numeric vector, the ids of the vertices or edges.  It is not
#' recycled, even if \code{value} is longer.
#' @param value Numeric vector, the new value(s) of the attributes, it will be
#' recycled if needed.
#' @return \code{get.graph.attribute}, \code{get.vertex.attribute} and
#' \code{get.edge.attribute} return an R object, or a list of R objects if
#' attributes of more vertices/edges are requested.
#' 
#' \code{set.graph.attribute}, \code{set.vertex.attribute},
#' \code{set.edge.attribute}, and also \code{remove.graph.attribute},
#' \code{remove.vertex.attribute} and \code{remove.edge.attribute} return a new
#' graph object with the updates/removes performed.
#' 
#' \code{list.graph.attributes}, \code{list.vertex.attributes} and
#' \code{list.edge.attributes} return a character vector, the names of the
#' attributes present.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{print.igraph}} can print attributes. See
#' \code{\link{attribute.combination}} for details on how igraph combines
#' attributes if several vertices or edges are mapped into one.
#' @keywords graphs
#' @examples
#' 
#' g <- graph.ring(10)
#' g <- set.graph.attribute(g, "name", "RING")
#' # It is the same as
#' g$name <- "RING"
#' g$name
#' 
#' g <- set.vertex.attribute(g, "color", value=c("red", "green"))
#' get.vertex.attribute(g, "color")
#' 
#' g <- set.edge.attribute(g, "weight", value=runif(ecount(g)))
#' get.edge.attribute(g, "weight")
#' 
#' # The following notation is more convenient
#' g <- graph.star(10)
#' 
#' V(g)$color <- c("red", "green")
#' V(g)$color
#' 
#' E(g)$weight <- runif(ecount(g))
#' E(g)$weight
#' 
#' str(g, g=TRUE, v=TRUE, e=TRUE)
#' 
#' # Setting all attributes at once
#' g2 <- graph.empty(10)
#' g2
#' 
#' graph.attributes(g2) <- graph.attributes(g)
#' vertex.attributes(g2) <- vertex.attributes(g)
#' edge.attributes(g2) <- list()
#' g2
#' 
#' edge.attributes(g2) <- list(weight=numeric())
#' g2
#' 
                            attributes=c("default", "all", "none")) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  if (!is.null(layers)) layers <- as.numeric(layers)-1
  hgap <- as.numeric(hgap)
  vgap <- as.numeric(vgap)
  maxiter <- as.integer(maxiter)
  if (is.null(weights) && "weight" %in% list.edge.attributes(graph)) { 
    weights <- E(graph)$weight 
  } 
  if (!is.null(weights) && any(!is.na(weights))) { 
    weights <- as.numeric(weights) 
  } else { 
    weights <- NULL 
  }
  attributes <- igraph.match.arg(attributes)
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_layout_sugiyama", graph, layers, hgap,
               vgap, maxiter, weights, PACKAGE="igraph")

  # Flip the y coordinates, more natural this way
  res$res[,2] <- max(res$res[,2]) - res$res[,2] + 1

  # Separate real and dummy vertices
  vc <- vcount(graph)
  res$layout <- res$res[seq_len(vc),]
  if (nrow(res$res)==vc) {
    res$layout.dummy <- matrix(nrow=0, ncol=2)
  } else {
    res$layout.dummy <- res$res[(vc+1):nrow(res$res),]
  }
  
  # Add some attributes to the extended graph
  E(res$extd_graph)$orig <- res$extd_to_orig_eids
  res$extd_to_orig_eids <- NULL

  res$extd_graph <- set.vertex.attribute(res$extd_graph, "dummy",
                                         value=c(rep(FALSE, vc),
                                           rep(TRUE, nrow(res$res)-vc)))

  res$extd_graph$layout <- rbind(res$layout, res$layout.dummy)

  if (attributes=="default" || attributes=="all") {
    if ("size" %in% list.vertex.attributes(graph)) {
      V(res$extd_graph)$size <- 0
      V(res$extd_graph)$size[ !V(res$extd_graph)$dummy ] <- V(graph)$size
    }
    if ("size2" %in% list.vertex.attributes(graph)) {
      V(res$extd_graph)$size2 <- 0
      V(res$extd_graph)$size2[ !V(res$extd_graph)$dummy ] <- V(graph)$size2
    }
    if ("shape" %in% list.vertex.attributes(graph)) {
      V(res$extd_graph)$shape <- "none"
      V(res$extd_graph)$shape[ !V(res$extd_graph)$dummy ] <- V(graph)$shape
    }
    if ("label" %in% list.vertex.attributes(graph)) {
      V(res$extd_graph)$label <- ""
      V(res$extd_graph)$label[ !V(res$extd_graph)$dummy ] <- V(graph)$label
    }
    if ("color" %in% list.vertex.attributes(graph)) {
      V(res$extd_graph)$color <- head(V(graph)$color, 1)
      V(res$extd_graph)$color[ !V(res$extd_graph)$dummy ] <- V(graph)$color
    }
    eetar <- get.edgelist(res$extd_graph, names=FALSE)[,2]
    E(res$extd_graph)$arrow.mode <- 0
    if ("arrow.mode" %in% list.edge.attributes(graph)) {
      E(res$extd_graph)$arrow.mode[ eetar <= vc ] <- E(graph)$arrow.mode
    } else {
      E(res$extd_graph)$arrow.mode[ eetar <= vc ] <- is.directed(graph) * 2
    }
    if ("arrow.size" %in% list.edge.attributes(graph)) {
      E(res$extd_graph)$arrow.size <- 0
      E(res$extd_graph)$arrow.size[ eetar <= vc ] <- E(graph)$arrow.size
    }
  }

  if (attributes=="all") {
    gatt <- setdiff(list.graph.attributes(graph), "layout")
    vatt <- setdiff(list.vertex.attributes(graph),
                    c("size", "size2", "shape", "label", "color"))
    eatt <- setdiff(list.edge.attributes(graph),
                    c("arrow.mode", "arrow.size"))
    for (ga in gatt) {
      res$extd_graph <- set.graph.attribute(res$extd_graph, ga,
                                            get.graph.attribute(graph, ga))
    }
    for (va in vatt) {
      notdummy <- which(!V(res$extd_graph)$dummy)
      res$extd_graph <- set.vertex.attribute(res$extd_graph, va,
                                             notdummy,
                                             get.vertex.attribute(graph, va))
    }
    for (ea in eatt) {
      eanew <- get.edge.attribute(graph, ea)[E(res$extd_graph)$orig]
      res$extd_graph <- set.edge.attribute(res$extd_graph, ea, value=eanew)
    }
  }
  
  res$res <- NULL
  res
}



#' Graph layout by multidimensional scaling
#' 
#' Multidimensional scaling of some distance matrix defined on the vertices of
#' a graph.
#' 
#' \code{layout.mds} uses metric multidimensional scaling for generating the
#' coordinates. Multidimensional scaling aims to place points from a higher
#' dimensional space in a (typically) 2 dimensional plane, so that the distance
#' between the points are kept as much as this is possible.
#' 
#' By default igraph uses the shortest path matrix as the distances between the
#' nodes, but the user can override this via the \code{dist} argument.
#' 
#' This function generates the layout separately for each graph component and
#' then merges them via \code{\link{layout.merge}}.
#' 
#' @param graph The input graph.
#' @param dist The distance matrix for the multidimensional scaling.  If
#' \code{NULL} (the default), then the unweighted shortest path matrix is used.
#' @param dim \code{layout.mds} supports dimensions up to the number of nodes
#' minus one, but only if the graph is connected; for unconnected graphs, the
#' only possible values is 2. This is because \code{layout.merge} only works in
#' 2D.
#' @param options This is currently ignored, as ARPACK is not used any more for
#' solving the eigenproblem
#' @return A numeric matrix with \code{dim} columns.
#' @author Tamas Nepusz \email{ntamas@@gmail.com} and Gabor Csardi
#' \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{layout}}, \code{\link{plot.igraph}}
#' @references Cox, T. F. and Cox, M. A. A. (2001) \emph{Multidimensional
#' Scaling}.  Second edition. Chapman and Hall.
#' @keywords graphs
#' @examples
#' 
#' g <- erdos.renyi.game(100, 2/100)
#' l <- layout.mds(g)
#' plot(g, layout=l, vertex.label=NA, vertex.size=3)
#' 
layout.mds <- function(graph, dist=NULL, dim=2,
                       options=igraph.arpack.default) {
  
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  if (!is.null(dist)) dist <- structure(as.double(dist), dim=dim(dist))
  dim <- as.integer(dim)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_layout_mds", graph, dist, dim,
        PACKAGE="igraph")

  res
}



#' The Fruchterman-Reingold layout algorithm
#' 
#' Place vertices on the plane using the force-directed layout algorithm by
#' Fruchterman and Reingold.
#' 
#' See the referenced paper below for the details of the algorithm.
#' 
#' This function was rewritten from scratch in igraph version 0.8.0.
#' 
#' @param graph The graph to lay out. Edge directions are ignored.
#' @param coords Optional starting positions for the vertices. If this argument
#' is not \code{NULL} then it should be an appropriate matrix of starting
#' coordinates.
#' @param dim Integer scalar, 2 or 3, the dimension of the layout.  Two
#' dimensional layouts are places on a plane, three dimensional ones in the 3d
#' space.
#' @param niter Integer scalar, the number of iterations to perform.
#' @param start.temp Real scalar, the start temperature. This is the maximum
#' amount of movement alloved along one axis, within one step, for a vertex.
#' Currently it is decreased linearly to zero during the iteration.
#' @param grid Character scalar, whether to use the faster, but less accurate
#' grid based implementation of the algorithm. By default (\dQuote{auto}), the
#' grid-based implementation is used if the graph has more than one thousand
#' vertices.
#' @param weights A vector giving edge weights. The \code{weight} edge
#' attribute is used by default, if present. If weights are given, then the
#' attraction along the edges will be multiplied by the given edge weights.
#' @param minx If not \code{NULL}, then it must be a numeric vector that gives
#' lower boundaries for the \sQuote{x} coordinates of the vertices. The length
#' of the vector must match the number of vertices in the graph.
#' @param maxx Similar to \code{minx}, but gives the upper boundaries.
#' @param miny Similar to \code{minx}, but gives the lower boundaries of the
#' \sQuote{y} coordinates.
#' @param maxy Similar to \code{minx}, but gives the upper boundaries of the
#' \sQuote{y} coordinates.
#' @param minz Similar to \code{minx}, but gives the lower boundaries of the
#' \sQuote{z} coordinates.
#' @param maxz Similar to \code{minx}, but gives the upper boundaries of the
#' \sQuote{z} coordinates.
#' @param coolexp,maxdelta,area,repulserad These arguments are not supported
#' from igraph version 0.8.0 and are ignored (with a warning).
#' @return A two- or three-column matrix, each row giving the coordinates of a
#' vertex, according to the ids of the vertex ids.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{layout.drl}}, \code{\link{layout.kamada.kawai}} for
#' other layout algorithms.
#' @references Fruchterman, T.M.J. and Reingold, E.M. (1991). Graph Drawing by
#' Force-directed Placement. \emph{Software - Practice and Experience},
#' 21(11):1129-1164.
#' @keywords graphs
#' @examples
#' 
#' # Fixing ego
#' g <- ba.game(20, m=2)
#' minC <- rep(-Inf, vcount(g))
#' maxC <- rep(Inf, vcount(g))
#' minC[1] <- maxC[1] <- 0
#' co <- layout.fruchterman.reingold(g, minx=minC, maxx=maxC,
#'                                   miny=minC, maxy=maxC)
#' co[1,]
#' plot(g, layout=co, vertex.size=30, edge.arrow.size=0.2,
#'      vertex.label=c("ego", rep("", vcount(g)-1)), rescale=FALSE,
#'      xlim=range(co[,1]), ylim=range(co[,2]), vertex.label.dist=0,
#'      vertex.label.color="red")
#' axis(1)
#' axis(2)
#' 
layout.fruchterman.reingold <- function(graph, coords=NULL, dim=2,
                            niter=500, start.temp=sqrt(vcount(graph)),
                            grid=c("auto", "grid", "nogrid"), weights=NULL,
                            minx=NULL, maxx=NULL, miny=NULL, maxy=NULL,
                            minz=NULL, maxz=NULL,
                            coolexp, maxdelta, area, repulserad) {

                                        # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  if (!is.null(coords)) {
    coords <- as.matrix(structure(as.double(coords), dim=dim(coords)))
  }
  dim <- as.integer(dim)
  if (dim != 2L && dim != 3L) {
    stop("Dimension must be two or three")
  }
  niter <- as.integer(niter)
  start.temp <- as.numeric(start.temp)

  grid <- igraph.match.arg(grid)
  grid <- switch(grid, "grid"=0L, "nogrid"=1L, "auto"=2L)

  if (is.null(weights) && "weight" %in% list.edge.attributes(graph)) {
    weights <- E(graph)$weight
  }
  if (!is.null(weights) && any(!is.na(weights))) {
    weights <- as.numeric(weights)
  } else {
    weights <- NULL
  }
  if (!is.null(minx)) minx <- as.numeric(minx)
  if (!is.null(maxx)) maxx <- as.numeric(maxx)
  if (!is.null(miny)) miny <- as.numeric(miny)
  if (!is.null(maxy)) maxy <- as.numeric(maxy)
  if (!is.null(minz)) minz <- as.numeric(minz)
  if (!is.null(maxz)) maxz <- as.numeric(maxz)
  if (!missing(coolexp)) {
    warning("Argument `coolexp' is deprecated and has no effect")
  }
  if (!missing(maxdelta)) {
    warning("Argument `maxdelta' is deprecated and has no effect")
  }
  if (!missing(area)) {
    warning("Argument `area' is deprecated and has no effect")
  }
  if (!missing(repulserad)) {
    warning("Argument `repulserad' is deprecated and has no effect")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  if (dim==2) {
    res <- .Call("R_igraph_layout_fruchterman_reingold", graph, coords,
                 niter, start.temp, weights, minx, maxx, miny, maxy, grid,
                 PACKAGE="igraph")
  } else {
    res <- .Call("R_igraph_layout_fruchterman_reingold_3d", graph, coords,
                 niter, start.temp, weights, minx, maxx, miny, maxy,
                 minz, maxz, PACKAGE="igraph")
  }
  res
}



#' The Kamada-Kawai layout algorithm
#' 
#' Place the vertices on the plane, or in the 3d space, based on a phyisical
#' model of springs.
#' 
#' See the referenced paper below for the details of the algorithm.
#' 
#' This function was rewritten from scratch in igraph version 0.8.0 and it
#' follows truthfully the original publication by Kamada and Kawai now.
#' 
#' @param graph The input graph. Edge directions are ignored.
#' @param coords If not \code{NULL}, then the starting coordinates should be
#' given here, in a two or three column matrix, depending on the \code{dim}
#' argument.
#' @param dim Integer scalar, 2 or 3, the dimension of the layout.  Two
#' dimensional layouts are places on a plane, three dimensional ones in the 3d
#' space.
#' @param maxiter The maximum number of iterations to perform. The algorithm
#' might terminate earlier, see the \code{epsilon} argument.
#' @param epsilon Numeric scalar, the algorithm terminates, if the maximal
#' delta is less than this. (See the reference below for what delta means.) If
#' you set this to zero, then the function always performs \code{maxiter}
#' iterations.
#' @param kkconst Numeric scalar, the Kamada-Kawai vertex attraction constant.
#' Typical (and default) value is the number of vertices.
#' @param weights Edge weights, larger values will result longer edges.
#' @param minx If not \code{NULL}, then it must be a numeric vector that gives
#' lower boundaries for the \sQuote{x} coordinates of the vertices. The length
#' of the vector must match the number of vertices in the graph.
#' @param maxx Similar to \code{minx}, but gives the upper boundaries.
#' @param miny Similar to \code{minx}, but gives the lower boundaries of the
#' \sQuote{y} coordinates.
#' @param maxy Similar to \code{minx}, but gives the upper boundaries of the
#' \sQuote{y} coordinates.
#' @param minz Similar to \code{minx}, but gives the lower boundaries of the
#' \sQuote{z} coordinates.
#' @param maxz Similar to \code{minx}, but gives the upper boundaries of the
#' \sQuote{z} coordinates.
#' @param niter,sigma,initemp,coolexp These arguments are not supported from
#' igraph version 0.8.0 and are ignored (with a warning).
#' @return A numeric matrix with two (dim=2) or three (dim=3) columns, and as
#' many rows as the number of vertices, the x, y and potentially z coordinates
#' of the vertices.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{layout.drl}}, \code{\link{plot.igraph}},
#' \code{\link{tkplot}}
#' @references Kamada, T. and Kawai, S.: An Algorithm for Drawing General
#' Undirected Graphs. \emph{Information Processing Letters}, 31/1, 7--15, 1989.
#' @keywords graphs
#' @examples
#' 
#' g <- graph.ring(10)
#' E(g)$weight <- rep(1:2, length.out=ecount(g))
#' plot(g, layout=layout.kamada.kawai, edge.label=E(g)$weight)
#' 
layout.kamada.kawai <- function(graph, coords=NULL, dim=2,
                                maxiter=50*vcount(graph),
                                epsilon=0.0, kkconst=vcount(graph),
                                weights=NULL, minx=NULL, maxx=NULL,
                                miny=NULL, maxy=NULL, minz=NULL, maxz=NULL,
                                niter, sigma, initemp, coolexp) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  if (!is.null(coords)) {
    coords <- as.matrix(structure(as.double(coords), dim=dim(coords)))
  }
  dim <- as.integer(dim)
  if (dim != 2L && dim != 3L) {
    stop("Dimension must be two or three")
  }

  maxiter <- as.integer(maxiter)
  epsilon <- as.numeric(epsilon)
  kkconst <- as.numeric(kkconst)
  if (is.null(weights) && "weight" %in% list.edge.attributes(graph)) { 
    weights <- E(graph)$weight 
  } 
  if (!is.null(weights) && any(!is.na(weights))) { 
    weights <- as.numeric(weights) 
  } else { 
    weights <- NULL 
  }
  if (!is.null(minx)) minx <- as.numeric(minx)
  if (!is.null(maxx)) maxx <- as.numeric(maxx)
  if (!is.null(miny)) miny <- as.numeric(miny)
  if (!is.null(maxy)) maxy <- as.numeric(maxy)
  if (!is.null(minz)) minz <- as.numeric(minz)
  if (!is.null(maxz)) maxz <- as.numeric(maxz)
  
  if (!missing(niter)) {
    warning("Argument `niter' is deprecated and has no effect")
  }
  if (!missing(sigma)) {
    warning("Argument `sigma' is deprecated and has no effect")
  }
  if (!missing(initemp)) {
    warning("Argument `initemp' is deprecated and has no effect")
  }
  if (!missing(coolexp)) {
    warning("Argument `coolexp' is deprecated and has no effect")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  if (dim == 2) {
    res <- .Call("R_igraph_layout_kamada_kawai", graph, coords, maxiter,
                 epsilon, kkconst, weights, minx, maxx, miny, maxy,
                 PACKAGE="igraph")
  } else {
    res <- .Call("R_igraph_layout_kamada_kawai_3d", graph, coords, maxiter,
                 epsilon, kkconst, weights, minx, maxx, miny, maxy, minz,
                 maxz, PACKAGE="igraph")
  }    

  res
}



#' The GEM layout algorithm
#' 
#' Place vertices on the plane using the GEM force-directed layout algorithm.
#' 
#' See the referenced paper below for the details of the algorithm.
#' 
#' @param graph The input graph. Edge directions are ignored.
#' @param coords If not \code{NULL}, then the starting coordinates should be
#' given here, in a two or three column matrix, depending on the \code{dim}
#' argument.
#' @param maxiter The maximum number of iterations to perform. Updating a
#' single vertex counts as an iteration.  A reasonable default is 40 * n * n,
#' where n is the number of vertices. The original paper suggests 4 * n * n,
#' but this usually only works if the other parameters are set up carefully.
#' @param temp.max The maximum allowed local temperature. A reasonable default
#' is the number of vertices.
#' @param temp.min The global temperature at which the algorithm terminates
#' (even before reaching \code{maxiter} iterations). A reasonable default is
#' 1/10.
#' @param temp.init Initial local temperature of all vertices. A reasonable
#' default is the square root of the number of vertices.
#' @return A numeric matrix with two columns, and as many rows as the number of
#' vertices.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{layout.fruchterman.reingold}},
#' \code{\link{plot.igraph}}, \code{\link{tkplot}}
#' @references Arne Frick, Andreas Ludwig, Heiko Mehldau: A Fast Adaptive
#' Layout Algorithm for Undirected Graphs, \emph{Proc. Graph Drawing 1994},
#' LNCS 894, pp. 388-403, 1995.
#' @keywords graphs
#' @examples
#' 
#' set.seed(42)
#' g <- graph.ring(10)
#' plot(g, layout=layout.gem)
#' 
layout.gem <- function(graph, coords=NULL, maxiter=40*vcount(graph)^2,
                       temp.max=vcount(graph), temp.min=1/10,
                       temp.init=sqrt(vcount(graph))) {
  
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  if (!is.null(coords)) {
    coords <- as.matrix(structure(as.double(coords), dim=dim(coords)))
    use.seed <- TRUE
  } else {
    coords <- matrix(ncol=2, nrow=0)
    use.seed <- FALSE
  }

  maxiter <- as.integer(maxiter)
  temp.max <- as.numeric(temp.max)
  temp.min <- as.numeric(temp.min)
  temp.init <- as.numeric(temp.init)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_layout_gem", graph, coords, use.seed, maxiter,
               temp.max, temp.min, temp.init,
               PACKAGE="igraph")

  res
}



#' The Davidson-Harel layout algorithm
#' 
#' Place vertices of a graph on the plane, according to the simulated annealing
#' algorithm by Davidson and Harel.
#' 
#' This function implements the algorithm by Davidson and Harel, see Ron
#' Davidson, David Harel: Drawing Graphs Nicely Using Simulated Annealing. ACM
#' Transactions on Graphics 15(4), pp. 301-331, 1996.
#' 
#' The algorithm uses simulated annealing and a sophisticated energy function,
#' which is unfortunately hard to parameterize for different graphs. The
#' original publication did not disclose any parameter values, and the ones
#' below were determined by experimentation.
#' 
#' The algorithm consists of two phases, an annealing phase, and a fine-tuning
#' phase. There is no simulated annealing in the second phase.
#' 
#' Our implementation tries to follow the original publication, as much as
#' possible. The only major difference is that coordinates are explicitly kept
#' within the bounds of the rectangle of the layout.
#' 
#' @param graph The graph to lay out. Edge directions are ignored.
#' @param coords Optional starting positions for the vertices. If this argument
#' is not \code{NULL} then it should be an appropriate matrix of starting
#' coordinates.
#' @param maxiter Number of iterations to perform in the first phase.
#' @param fineiter Number of iterations in the fine tuning phase.
#' @param cool.fact Cooling factor.
#' @param weight.node.dist Weight for the node-node distances component of the
#' energy function.
#' @param weight.border Weight for the distance from the border component of
#' the energy function. It can be set to zero, if vertices are allowed to sit
#' on the border.
#' @param weight.edge.lengths Weight for the edge length component of the
#' energy function.
#' @param weight.edge.crossings Weight for the edge crossing component of the
#' energy function.
#' @param weight.node.edge.dist Weight for the node-edge distance component of
#' the energy function.
#' @return A two- or three-column matrix, each row giving the coordinates of a
#' vertex, according to the ids of the vertex ids.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{layout.fruchterman.reingold}},
#' \code{\link{layout.kamada.kawai}} for other layout algorithms.
#' @references Ron Davidson, David Harel: Drawing Graphs Nicely Using Simulated
#' Annealing. \emph{ACM Transactions on Graphics} 15(4), pp. 301-331, 1996.
#' @examples
#' 
#' set.seed(42)
#' L <- layout.davidson.harel
#' ## Figures from the paper
#' g_1b <- graph.star(19, mode="undirected") + path(c(2:19, 2)) +
#'   path(c(seq(2, 18, by=2), 2))
#' plot(g_1b, layout=L)
#' 
#' g_2 <- graph.lattice(c(8, 3)) + edges(1,8, 9,16, 17,24)
#' plot(g_2, layout=L)
#' 
#' g_3 <- graph.empty(n=70)
#' plot(g_3, layout=L)
#' 
#' g_4 <- graph.empty(n=70, directed=FALSE) + edges(1:70)
#' plot(g_4, layout=L, vertex.size=5, vertex.label=NA)
#' 
#' g_5a <- graph.ring(24)
#' plot(g_5a, layout=L, vertex.size=5, vertex.label=NA)
#' 
#' g_5b <- graph.ring(40)
#' plot(g_5b, layout=L, vertex.size=5, vertex.label=NA)
#' 
#' g_6 <- graph.lattice(c(2,2,2))
#' plot(g_6, layout=L)
#' 
#' g_7 <- graph.formula(1:3:5 -- 2:4:6)
#' plot(g_7, layout=L, vertex.label=V(g_7)$name)
#' 
#' g_8 <- graph.ring(5) + graph.ring(10) + graph.ring(5) +
#'   edges(1,6, 2,8, 3, 10, 4,12, 5,14,
#'         7,16, 9,17, 11,18, 13,19, 15,20)
#' plot(g_8, layout=L, vertex.size=5, vertex.label=NA)
#' 
#' g_9 <- graph.lattice(c(3,2,2))
#' plot(g_9, layout=L, vertex.size=5, vertex.label=NA)
#' 
#' g_10 <- graph.lattice(c(6,6))
#' plot(g_10, layout=L, vertex.size=5, vertex.label=NA)
#' 
#' g_11a <- graph.tree(31, 2, mode="undirected")
#' plot(g_11a, layout=L, vertex.size=5, vertex.label=NA)
#' 
#' g_11b <- graph.tree(21, 4, mode="undirected")
#' plot(g_11b, layout=L, vertex.size=5, vertex.label=NA)
#' 
#' g_12 <- graph.empty(n=37, directed=FALSE) +
#'   path(1:5,10,22,31,37:33,27,16,6,1) + path(6,7,11,9,10) + path(16:22) +
#'   path(27:31) + path(2,7,18,28,34) + path(3,8,11,19,29,32,35) +
#'   path(4,9,20,30,36) + path(1,7,12,14,19,24,26,30,37) +
#'   path(5,9,13,15,19,23,25,28,33) + path(3,12,16,25,35,26,22,13,3)
#' plot(g_12,  layout=L, vertex.size=5, vertex.label=NA)
#' 
layout.davidson.harel <- function(graph, coords=NULL, maxiter=10,
           fineiter=max(10, log2(vcount(graph))), cool.fact=0.75,
           weight.node.dist=1.0, weight.border=0.0,
           weight.edge.lengths=graph.density(graph) / 10,
           weight.edge.crossings=1.0 - sqrt(graph.density(graph)),
           weight.node.edge.dist=0.2 * (1-graph.density(graph))) {
  
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  if (!is.null(coords)) {
    coords <- as.matrix(structure(as.double(coords), dim=dim(coords)))
    use.seed <- TRUE
  } else {
    coords <- matrix(ncol=2, nrow=0)
    use.seed <- FALSE
  }
  maxiter <- as.integer(maxiter)
  fineiter <- as.integer(fineiter)
  cool.fact <- as.numeric(cool.fact)
  weight.node.dist <- as.numeric(weight.node.dist)
  weight.border <- as.numeric(weight.border)
  weight.edge.lengths <- as.numeric(weight.edge.lengths)
  weight.edge.crossings <- as.numeric(weight.edge.crossings)
  weight.node.edge.dist <- as.numeric(weight.node.edge.dist)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_layout_davidson_harel", graph, coords, use.seed,
               maxiter, fineiter, cool.fact, weight.node.dist,
               weight.border, weight.edge.lengths, weight.edge.crossings,
               weight.node.edge.dist, PACKAGE="igraph")

  res
}

