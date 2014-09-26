#   IGraph R package
#   Copyright (C) 2005-2012  Gabor Csardi <csardi.gabor@gmail.com>
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
# Connected components, subgraphs, kinda
###################################################################

no.clusters <- function(graph, mode=c("weak", "strong")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "weak"=1, "strong"=2)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_no_clusters", graph, as.numeric(mode),
        PACKAGE="igraph")
}

#' @rdname clusters
#' @param cumulative Logical, if TRUE the cumulative distirubution (relative
#' frequency) is calculated.
#' @param mul.size Logical. If TRUE the relative frequencies will be multiplied
#' by the cluster sizes.

cluster.distribution <- function(graph, cumulative=FALSE, mul.size=FALSE,
                                 ...) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  
  cs <- clusters(graph, ...)$csize;
  hi <- hist(cs, -1:max(cs), plot=FALSE)$density
  if (mul.size) {
    hi <- hi*1:max(cs)
    hi <- hi/sum(hi)
  }
  if (!cumulative) {
    res <- hi
  } else {
    res <- rev(cumsum(rev(hi)));
  }
  
  res
}

is.connected <- function(graph, mode=c("weak", "strong")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "weak"=1, "strong"=2)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_is_connected", graph, as.numeric(mode),
        PACKAGE="igraph")
}



#' Decompose a graph into components
#' 
#' Creates a separate graph for each component of a graph.
#' 
#' 
#' @param graph The original graph.
#' @param mode Character constant giving the type of the components, wither
#' \code{weak} for weakly connected components or \code{strong} for strongly
#' connected components.
#' @param max.comps The maximum number of components to return. The first
#' \code{max.comps} components will be returned (which hold at least
#' \code{min.vertices} vertices, see the next parameter), the others will be
#' ignored. Supply \code{NA} here if you don't want to limit the number of
#' components.
#' @param min.vertices The minimum number of vertices a component should
#' contain in order to place it in the result list. Eg. supply 2 here to ignore
#' isolate vertices.
#' @return A list of graph objects.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{is.connected}} to decide whether a graph is connected,
#' \code{\link{clusters}} to calculate the connected components of a graph.
#' @keywords graphs
#' @examples
#' 
#' # the diameter of each component in a random graph
#' g <- erdos.renyi.game(1000, 1/1000)
#' comps <- decompose.graph(g, min.vertices=2)
#' sapply(comps, diameter)
#' 
decompose.graph <- function(graph, mode=c("weak", "strong"), max.comps=NA,
                      min.vertices=0) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "weak"=1, "strong"=2)

  if (is.na(max.comps)) {
    max.comps=-1
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_decompose", graph, as.numeric(mode),
        as.numeric(max.comps), as.numeric(min.vertices),
        PACKAGE="igraph"
        )
}
