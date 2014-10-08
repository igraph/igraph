#   IGraph R package
#   Copyright (C) 2009-2012  Gabor Csardi <csardi.gabor@gmail.com>
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



#' Project a bipartite graph
#' 
#' A bipartite graph is projected into two one-mode networks
#' 
#' Bipartite graphs have a \code{type} vertex attribute in igraph, this is
#' boolean and \code{FALSE} for the vertices of the first kind and \code{TRUE}
#' for vertices of the second kind.
#' 
#' \code{bipartite_projection_size} calculates the number of vertices and edges
#' in the two projections of the bipartite graphs, without calculating the
#' projections themselves. This is useful to check how much memory the
#' projections would need if you have a large bipartite graph.
#' 
#' \code{bipartite_projection} calculates the actual projections.  You can use
#' the \code{probe1} argument to specify the order of the projections in the
#' result. By default vertex type \code{FALSE} is the first and \code{TRUE} is
#' the second.
#' 
#' \code{bipartite_projection} keeps vertex attributes.
#' 
#' @aliases bipartite.projection bipartite.projection.size bipartite_projection_size bipartite_projection
#' @param graph The input graph. It can be directed, but edge directions are
#' ignored during the computation.
#' @param types An optional vertex type vector to use instead of the
#' \sQuote{\code{type}} vertex attribute. You must supply this argument if the
#' graph has no \sQuote{\code{type}} vertex attribute.
#' @param multiplicity If \code{TRUE}, then igraph keeps the multiplicity of
#' the edges as an edge attribute. E.g. if there is an A-C-B and also an A-D-B
#' triple in the bipartite graph (but no more X, such that A-X-B is also in the
#' graph), then the multiplicity of the A-B edge in the projection will be 2.
#' @param probe1 This argument can be used to specify the order of the
#' projections in the resulting list. If given, then it is considered as a
#' vertex id (or a symbolic vertex name); the projection containing this vertex
#' will be the first one in the result list.  This argument is ignored if only
#' one projection is requested in argument \code{which}.
#' @param which A character scalar to specify which projection(s) to calculate.
#' The default is to calculate both.
#' @param remove.type Logical scalar, whether to remove the \code{type} vertex
#' attribute from the projections. This makes sense because these graphs are
#' not bipartite any more. However if you want to combine them with each other
#' (or other bipartite graphs), then it is worth keeping this attribute. By
#' default it will be removed.
#' @return A list of two undirected graphs. See details above.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @export
#' @keywords graphs
#' @examples
#' 
#' ## Projection of a full bipartite graph is a full graph
#' g <- full_bipartite_graph(10,5)
#' proj <- bipartite_projection(g)
#' graph.isomorphic(proj[[1]], full_graph(10))
#' graph.isomorphic(proj[[2]], full_graph(5))
#' 
#' ## The projection keeps the vertex attributes
#' M <- matrix(0, nr=5, nc=3)
#' rownames(M) <- c("Alice", "Bob", "Cecil", "Dan", "Ethel")
#' colnames(M) <- c("Party", "Skiing", "Badminton")
#' M[] <- sample(0:1, length(M), replace=TRUE)
#' M
#' g2 <- graph_from_incidence_matrix(M)
#' g2$name <- "Event network"
#' proj2 <- bipartite_projection(g2)
#' print(proj2[[1]], g=TRUE, e=TRUE)
#' print(proj2[[2]], g=TRUE, e=TRUE)
#' 
bipartite_projection <- function(graph, types=NULL,
                                 multiplicity=TRUE, probe1=NULL,
				 which=c("both", "true", "false"),
                                 remove.type=TRUE) {
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  if (is.null(types) && "type" %in% vertex_attr_names(graph)) { 
  types <- V(graph)$type 
  } 
  if (!is.null(types)) {
    if (!is.logical(types)) {
      warning("vertex types converted to logical")
    }
    types <- as.logical(types)
    if (any(is.na(types))) {
      stop("`NA' is not allowed in vertex types")
    }
  } else { 
  stop("Not a bipartite graph, supply `types' argument") 
  }
  if (!is.null(probe1)) {
    probe1 <- as.igraph.vs(graph, probe1)-1
  } else {
    probe1 <- -1
  }
  which <- switch(igraph.match.arg(which), "both"=0L, "false"=1L,
                  "true"=2L)
  if (which != "both" && probe1 != -1) {
    warning("`probe1' ignored if only one projection is requested")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_bipartite_projection", graph, types,
               as.integer(probe1), which, PACKAGE="igraph")
  if (remove.type) {
    if (is_igraph(res[[1]])) {
      res[[1]] <- delete_vertex_attr(res[[1]], "type")
    }
    if (is_igraph(res[[2]])) {
      res[[2]] <- delete_vertex_attr(res[[2]], "type")
    }
  }
  if (which == 0L) {
    if (multiplicity) {
      E(res[[1]])$weight <- res[[3]]
      E(res[[2]])$weight <- res[[4]]
    }
    res[1:2]
  } else if (which == 1L) {
    if (multiplicity) { E(res[[1]])$weight <- res[[3]] }
    res[[1]]
  } else {
    if (multiplicity) { E(res[[2]])$weight <- res[[4]] }
    res[[2]]
  }
}
