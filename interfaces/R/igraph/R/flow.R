#   IGraph R package
#   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

min_cut <- function(graph, source=NULL, target=NULL, capacity=NULL,
                         value.only=TRUE) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  if (is.null(capacity)) {
    if ("capacity" %in% edge_attr_names(graph)) {
      capacity <- E(graph)$capacity
    }
  }
  if (is.null(source) && !is.null(target) ||
      is.null(target) && !is.null(source)) {
    stop("Please give both source and target or neither")
  }
  if (!is.null(capacity)) {
    capacity <- as.numeric(capacity)
  }

  value.only <- as.logical(value.only)
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  if (is.null(target) && is.null(source)) {
    if (value.only) {
      res <- .Call("R_igraph_mincut_value", graph, capacity,
                   PACKAGE="igraph")
    } else {
      res <- .Call("R_igraph_mincut", graph, capacity,
                   PACKAGE="igraph")
      res$cut <- res$cut + 1
      res$partition1 <- res$partition1 + 1
      res$partition2 <- res$partition2 + 1
      res
    }
  } else {
    if (value.only) {
      res <- .Call("R_igraph_st_mincut_value", graph,
                   as.igraph.vs(graph, source)-1,
                   as.igraph.vs(graph, target)-1, capacity,
                   PACKAGE="igraph")
    } else {
      stop("Calculating minimum s-t cuts is not implemented yet")
    }
  }
  res
}



#' Vertex connectivity.
#' 
#' The vertex connectivity of a graph or two vertices, this is recently also
#' called group cohesion.
#' 
#' The vertex connectivity of two vertices (\code{source} and \code{target}) in
#' a directed graph is the minimum number of vertices needed to remove from the
#' graph to eliminate all (directed) paths from \code{source} to \code{target}.
#' \code{vertex_connectivity} calculates this quantity if both the
#' \code{source} and \code{target} arguments are given and they're not
#' \code{NULL}.
#' 
#' The vertex connectivity of a graph is the minimum vertex connectivity of all
#' (ordered) pairs of vertices in the graph. In other words this is the minimum
#' number of vertices needed to remove to make the graph not strongly
#' connected. (If the graph is not strongly connected then this is zero.)
#' \code{vertex_connectivity} calculates this quantitty if neither the
#' \code{source} nor \code{target} arguments are given. (Ie. they are both
#' \code{NULL}.)
#' 
#' A set of vertex disjoint directed paths from \code{source} to \code{vertex}
#' is a set of directed paths between them whose vertices do not contain common
#' vertices (apart from \code{source} and \code{target}). The maximum number of
#' vertex disjoint paths between two vertices is the same as their vertex
#' connectivity in most cases (if the two vertices are not connected by an
#' edge).
#' 
#' The cohesion of a graph (as defined by White and Harary, see references), is
#' the vertex connectivity of the graph. This is calculated by
#' \code{cohesion}.
#' 
#' These three functions essentially calculate the same measure(s), more
#' precisely \code{vertex_connectivity} is the most general, the other two are
#' included only for the ease of using more descriptive function names.
#' 
#' @aliases vertex.connectivity vertex.disjoint.paths cohesion vertex_connectivity
#'   vertex_disjoint_paths graph.cohesion
#' @param graph,x The input graph.
#' @param source The id of the source vertex, for \code{vertex_connectivity} it
#' can be \code{NULL}, see details below.
#' @param target The id of the target vertex, for \code{vertex_connectivity} it
#' can be \code{NULL}, see details below.
#' @param checks Logical constant. Whether to check that the graph is connected
#' and also the degree of the vertices. If the graph is not (strongly)
#' connected then the connectivity is obviously zero. Otherwise if the minimum
#' degree is one then the vertex connectivity is also one. It is a good idea to
#' perform these checks, as they can be done quickly compared to the
#' connectivity calculation itself.  They were suggested by Peter McMahan,
#' thanks Peter.
#' @param ... Ignored.
#' @return A scalar real value.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{max_flow}}, \code{\link{edge_connectivity}},
#' \code{\link{edge_disjoint_paths}}, \code{\link{adhesion}}
#' @references White, Douglas R and Frank Harary 2001. The Cohesiveness of
#' Blocks In Social Networks: Node Connectivity and Conditional Density.
#' \emph{Sociological Methodology} 31 (1) : 305-359.
#' @keywords graphs
#' @examples
#' 
#' g <- barabasi.game(100, m=1)
#' g <- delete_edges(g, E(g)[ 100 %--% 1 ])
#' g2 <- barabasi.game(100, m=5)
#' g2 <- delete_edges(g2, E(g2)[ 100 %--% 1])
#' vertex_connectivity(g, 100, 1)
#' vertex_connectivity(g2, 100, 1)
#' vertex_disjoint_paths(g2, 100, 1)
#' 
#' g <- sample_gnp(50, 5/50)
#' g <- as.directed(g)
#' g <- induced_subgraph(g, subcomponent(g, 1))
#' cohesion(g)
#' 
vertex_connectivity <- function(graph, source=NULL, target=NULL, checks=TRUE) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  if (is.null(source) && is.null(target)) {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    .Call("R_igraph_vertex_connectivity", graph, as.logical(checks),
          PACKAGE="igraph")
  } else if (!is.null(source) && !is.null(target)) {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    .Call("R_igraph_st_vertex_connectivity", graph, as.igraph.vs(graph, source)-1,
          as.igraph.vs(graph, target)-1,
          PACKAGE="igraph")
  } else {
    stop("either give both source and target or neither")
  }
}



#' Edge connectivity.
#' 
#' The edge connectivity of a graph or two vertices, this is recently also
#' called group adhesion.
#' 
#' The edge connectivity of a pair of vertices (\code{source} and
#' \code{target}) is the minimum number of edges needed to remove to eliminate
#' all (directed) paths from \code{source} to \code{target}.
#' \code{edge_connectivity} calculates this quantity if both the \code{source}
#' and \code{target} arguments are given (and not \code{NULL}).
#' 
#' The edge connectivity of a graph is the minimum of the edge connectivity of
#' every (ordered) pair of vertices in the graph.  \code{edge_connectivity}
#' calculates this quantity if neither the \code{source} nor the \code{target}
#' arguments are given (ie. they are both \code{NULL}).
#' 
#' A set of edge disjoint paths between two vertices is a set of paths between
#' them containing no common edges. The maximum number of edge disjoint paths
#' between two vertices is the same as their edge connectivity.
#' 
#' The adhesion of a graph is the minimum number of edges needed to remove to
#' obtain a graph which is not strongly connected. This is the same as the edge
#' connectivity of the graph.
#' 
#' The three functions documented on this page calculate similar properties,
#' more precisely the most general is \code{edge_connectivity}, the others are
#' included only for having more descriptive function names.
#' 
#' @aliases edge.connectivity edge_disjoint_paths graph.adhesion adhesion
#'   edge_connectivity edge.disjoint.paths
#' @param graph The input graph.
#' @param source The id of the source vertex, for \code{edge_connectivity} it
#' can be \code{NULL}, see details below.
#' @param target The id of the target vertex, for \code{edge_connectivity} it
#' can be \code{NULL}, see details below.
#' @param checks Logical constant. Whether to check that the graph is connected
#' and also the degree of the vertices. If the graph is not (strongly)
#' connected then the connectivity is obviously zero. Otherwise if the minimum
#' degree is one then the edge connectivity is also one. It is a good idea to
#' perform these checks, as they can be done quickly compared to the
#' connectivity calculation itself.  They were suggested by Peter McMahan,
#' thanks Peter.
#' @return A scalar real value.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{max_flow}}, \code{\link{vertex_connectivity}},
#' \code{\link{vertex_disjoint_paths}}, \code{\link{cohesion}}
#' @references Douglas R. White and Frank Harary: The cohesiveness of blocks in
#' social networks: node connectivity and conditional density, TODO: citation
#' @keywords graphs
#' @examples
#' 
#' g <- barabasi.game(100, m=1)
#' g2 <- barabasi.game(100, m=5)
#' edge_connectivity(g, 100, 1)
#' edge_connectivity(g2, 100, 1)
#' edge_disjoint_paths(g2, 100, 1)
#' 
#' g <- sample_gnp(50, 5/50)
#' g <- as.directed(g)
#' g <- induced_subgraph(g, subcomponent(g, 1))
#' adhesion(g)
#' 
edge_connectivity <- function(graph, source=NULL, target=NULL, checks=TRUE) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  if (is.null(source) && is.null(target)) {    
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    .Call("R_igraph_edge_connectivity", graph, as.logical(checks),
          PACKAGE="igraph")
  } else if (!is.null(source) && !is.null(target)) {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    .Call("R_igraph_st_edge_connectivity", graph,
          as.igraph.vs(graph, source)-1, as.igraph.vs(graph, target)-1,
          PACKAGE="igraph")
  } else {
    stop("either give both source and target or neither")
  }
}

edge_disjoint_paths <- function(graph, source, target) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_edge_disjoint_paths", graph,
        as.igraph.vs(graph, source)-1, as.igraph.vs(graph, target)-1,
        PACKAGE="igraph")
}

vertex_disjoint_paths <- function(graph, source=NULL, target=NULL) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_vertex_disjoint_paths", graph, as.igraph.vs(graph, source)-1,
        as.igraph.vs(graph, target)-1,
        PACKAGE="igraph")
}

adhesion <- function(graph, checks=TRUE) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_adhesion", graph, as.logical(checks),
        PACKAGE="igraph")
}

#' @rdname vertex_connectivity
#' @method cohesion igraph

cohesion.igraph <- function(x, checks=TRUE, ...) {

  if (!is_igraph(x)) {
    stop("Not a graph object")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_cohesion", x, as.logical(checks),
        PACKAGE="igraph")
}

