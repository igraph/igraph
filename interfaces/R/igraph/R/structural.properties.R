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
# Structural properties
###################################################################



#' Diameter of a graph
#' 
#' The diameter of a graph is the length of the longest geodesic.
#' 
#' The diameter is calculated by using a breadth-first search like method.
#' 
#' \code{get_diameter} returns a path with the actual diameter. If there are
#' many shortest paths of the length of the diameter, then it returns the first
#' one found.
#' 
#' \code{farthest.points} returns two vertex ids, the vertices which are
#' connected by the diameter path.
#' 
#' @aliases diameter get.diameter farthest.nodes farthest_vertices get_diameter
#' @param graph The graph to analyze.
#' @param directed Logical, whether directed or undirected paths are to be
#' considered. This is ignored for undirected graphs.
#' @param unconnected Logical, what to do if the graph is unconnected. If
#' FALSE, the function will return a number that is one larger the largest
#' possible diameter, which is always the number of vertices. If TRUE, the
#' diameters of the connected components will be calculated and the largest one
#' will be returned.
#' @param weights Optional positive weight vector for calculating weighted
#' distances. If the graph has a \code{weight} edge attribute, then this is
#' used by default.
#' @return A numeric constant for \code{diameter}, a numeric vector for
#' \code{get_diameter} and a numeric vector of length two for
#' \code{farthest_vertices}.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{distances}}
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- make_ring(10)
#' g2 <- delete_edges(g, c(1,2,1,10))
#' diameter(g2, unconnected=TRUE)
#' diameter(g2, unconnected=FALSE)
#' 
#' ## Weighted diameter
#' set.seed(1)
#' g <- make_ring(10)
#' E(g)$weight <- sample(seq_len(ecount(g)))
#' diameter(g)
#' get_diameter(g)
#' diameter(g, weights=NA)
#' get_diameter(g, weights=NA)
#' 
diameter <- function(graph, directed=TRUE, unconnected=TRUE, weights=NULL) {
  
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  if (is.null(weights) && "weight" %in% edge_attr_names(graph)) {
    weights <- E(graph)$weight
  }
  if (!is.null(weights) && any(!is.na(weights))) {
    weights <- as.numeric(weights)
  } else {
    weights <- NULL
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_diameter", graph, as.logical(directed),
        as.logical(unconnected), weights,
        PACKAGE="igraph")
}

#' @export

get_diameter <- function(graph, directed=TRUE, unconnected=TRUE,
                         weights=NULL) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  if (is.null(weights) && "weight" %in% edge_attr_names(graph)) {
    weights <- E(graph)$weight
  }
  if (!is.null(weights) && any(!is.na(weights))) {
    weights <- as.numeric(weights)
  } else {
    weights <- NULL
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_get_diameter", graph, as.logical(directed),
               as.logical(unconnected), weights,
               PACKAGE="igraph")
  res + 1
}

#' @export

farthest_vertices <- function(graph, directed=TRUE, unconnected=TRUE,
                           weights=NULL) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  if (is.null(weights) && "weight" %in% edge_attr_names(graph)) {
    weights <- E(graph)$weight
  }
  if (!is.null(weights) && any(!is.na(weights))) {
    weights <- as.numeric(weights)
  } else {
    weights <- NULL
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_farthest_points", graph, as.logical(directed),
               as.logical(unconnected), weights,
               PACKAGE="igraph")
  res[1:2] <- res[1:2] + 1
  res
}       

#' @export

mean_distance <- function(graph, directed=TRUE, unconnected=TRUE) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_average_path_length", graph, as.logical(directed),
        as.logical(unconnected),
        PACKAGE="igraph")
}



#' Degree and degree distribution of the vertices
#' 
#' The degree of a vertex is its most basic structural property, the number of
#' its adjacent edges.
#' 
#' 
#' @aliases degree degree.distribution degree_distribution
#' @param graph The graph to analyze.
#' @param v The ids of vertices of which the degree will be calculated.
#' @param mode Character string, \dQuote{out} for out-degree, \dQuote{in} for
#' in-degree or \dQuote{total} for the sum of the two. For undirected graphs
#' this argument is ignored. \dQuote{all} is a synonym of \dQuote{total}.
#' @param loops Logical; whether the loop edges are also counted.
#' @param normalized Logical scalar, whether to normalize the degree.  If
#' \code{TRUE} then the result is divided by \eqn{n-1}, where \eqn{n} is the
#' number of vertices in the graph.
#' @param \dots Additional arguments to pass to \code{degree}, eg. \code{mode}
#' is useful but also \code{v} and \code{loops} make sense.
#' @return For \code{degree} a numeric vector of the same length as argument
#' \code{v}.
#' 
#' For \code{degree_distribution} a numeric vector of the same length as the
#' maximum degree plus one. The first element is the relative frequency zero
#' degree vertices, the second vertices with degree one, etc.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @keywords graphs
#' @export
#' @examples
#' 
#' g <- make_ring(10)
#' degree(g)
#' g2 <- sample_gnp(1000, 10/1000)
#' degree_distribution(g2)
#' 
degree <- function(graph, v=V(graph),
                   mode=c("all", "out", "in", "total"), loops=TRUE,
                   normalized=FALSE){
  
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  v <- as.igraph.vs(graph, v)
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3, "total"=3)
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_degree", graph, v-1,
               as.numeric(mode), as.logical(loops), PACKAGE="igraph")
  if (normalized) { res <- res / (vcount(graph)-1) }
  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    names(res) <- V(graph)$name[v]
  }
  res
}

#' @rdname degree
#' @param cumulative Logical; whether the cumulative degree distribution is to
#' be calculated.
#' @export

degree_distribution <- function(graph, cumulative=FALSE, ...) {
  
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  cs <- degree(graph, ...)
  hi <- hist(cs, -1:max(cs), plot=FALSE)$density
  if (!cumulative) {
    res <- hi
  } else {
    res <- rev(cumsum(rev(hi)))
  }
  
  res
}



#' Shortest (directed or undirected) paths between vertices
#' 
#' \code{distances} calculates the length of all the shortest paths from
#' or to the vertices in the network. \code{shortest_paths} calculates one
#' shortest path (the path itself, and not just its length) from or to the
#' given vertex.
#' 
#' The shortest path, or geodesic between two pair of vertices is a path with
#' the minimal number of vertices. The functions documented in this manual page
#' all calculate shortest paths between vertex pairs.
#' 
#' \code{distances} calculates the lengths of pairwise shortest paths from
#' a set of vertices (\code{from}) to another set of vertices (\code{to}). It
#' uses different algorithms, depending on the \code{argorithm} argument and
#' the \code{weight} edge attribute of the graph. The implemented algorithms
#' are breadth-first search (\sQuote{\code{unweighted}}), this only works for
#' unweighted graphs; the Dijkstra algorithm (\sQuote{\code{dijkstra}}), this
#' works for graphs with non-negative edge weights; the Bellman-Ford algorithm
#' (\sQuote{\code{bellman-ford}}), and Johnson's algorithm
#' (\sQuote{\code{"johnson"}}). The latter two algorithms work with arbitrary
#' edge weights, but (naturally) only for graphs that don't have a negative
#' cycle.
#' 
#' igraph can choose automatically between algorithms, and chooses the most
#' efficient one that is appropriate for the supplied weights (if any). For
#' automatic algorithm selection, supply \sQuote{\code{automatic}} as the
#' \code{algorithm} argument. (This is also the default.)
#' 
#' \code{shortest_paths} calculates a single shortest path (i.e. the path
#' itself, not just its length) between the source vertex given in \code{from},
#' to the target vertices given in \code{to}. \code{shortest_paths} uses
#' breadth-first search for unweighted graphs and Dijkstra's algorithm for
#' weighted graphs. The latter only works if the edge weights are non-negative.
#' 
#' \code{all_shortest_paths} calculates \emph{all} shortest paths between
#' pairs of vertices. More precisely, between the \code{from} vertex to the
#' vertices given in \code{to}. It uses a breadth-first search for unweighted
#' graphs and Dijkstra's algorithm for weighted ones. The latter only supports
#' non-negative edge weights.
#' 
#' \code{mean_distance} calculates the average path length in a graph, by
#' calculating the shortest paths between all pairs of vertices (both ways for
#' directed graphs). This function does not consider edge weights currently and
#' uses a breadth-first search.
#' 
#' \code{distance_table} calculates a histogram, by calculating the shortest
#' path length between each pair of vertices. For directed graphs both
#' directions are considered, so every pair of vertices appears twice in the
#' histogram.
#' 
#' @aliases shortest.paths get.shortest.paths get.all.shortest.paths distances
#' mean_distance distance_table average.path.length path.length.hist
#' all_shortest_paths shortest_paths
#' @param graph The graph to work on.
#' @param v Numeric vector, the vertices from which the shortest paths will be
#' calculated.
#' @param to Numeric vector, the vertices to which the shortest paths will be
#' calculated. By default it includes all vertices. Note that for
#' \code{distances} every vertex must be included here at most once. (This
#' is not required for \code{shortest_paths}.
#' @param mode Character constant, gives whether the shortest paths to or from
#' the given vertices should be calculated for directed graphs. If \code{out}
#' then the shortest paths \emph{from} the vertex, if \code{in} then \emph{to}
#' it will be considered. If \code{all}, the default, then the corresponding
#' undirected graph will be used, ie. not directed paths are searched. This
#' argument is ignored for undirected graphs.
#' @param weights Possibly a numeric vector giving edge weights. If this is
#' \code{NULL} and the graph has a \code{weight} edge attribute, then the
#' attribute is used. If this is \code{NA} then no weights are used (even if
#' the graph has a \code{weight} attribute).
#' @param algorithm Which algorithm to use for the calculation. By default
#' igraph tries to select the fastest suitable algorithm. If there are no
#' weights, then an unweighted breadth-first search is used, otherwise if all
#' weights are positive, then Dijkstra's algorithm is used. If there are
#' negative weights and we do the calculation for more than 100 sources, then
#' Johnson's algorithm is used. Otherwise the Bellman-Ford algorithm is used.
#' You can override igraph's choice by explicitly giving this parameter. Note
#' that the igraph C core might still override your choice in obvious cases,
#' i.e. if there are no edge weights, then the unweighted algorithm will be
#' used, regardless of this argument.
#' @return For \code{distances} a numeric matrix with \code{length(to)}
#' columns and \code{length(v)} rows. The shortest path length from a vertex to
#' itself is always zero. For unreachable vertices \code{Inf} is included.
#' 
#' For \code{shortest_paths} a named list with four entries is returned:
#' \item{vpath}{This itself is a list, of length \code{length(to)}; list
#' element \code{i} contains the vertex ids on the path from vertex \code{from}
#' to vertex \code{to[i]} (or the other way for directed graphs depending on
#' the \code{mode} argument). The vector also contains \code{from} and \code{i}
#' as the first and last elements. If \code{from} is the same as \code{i} then
#' it is only included once. If there is no path between two vertices then a
#' numeric vector of length zero is returned as the list element. If this
#' output is not requested in the \code{output} argument, then it will be
#' \code{NULL}.} \item{epath}{This is a list similar to \code{vpath}, but the
#' vectors of the list contain the edge ids along the shortest paths, instead
#' of the vertex ids. This entry is set to \code{NULL} if it is not requested
#' in the \code{output} argument.} \item{predecessors}{Numeric vector, the
#' predecessor of each vertex in the \code{to} argument, or \code{NULL} if it
#' was not requested.} \item{inbound_edges}{Numeric vector, the inbound edge
#' for each vertex, or \code{NULL}, if it was not requested.}
#' 
#' For \code{all_shortest_paths} a list is returned, each list element
#' contains a shortest path from \code{from} to a vertex in \code{to}. The
#' shortest paths to the same vertex are collected into consecutive elements of
#' the list.
#' 
#' For \code{mean_distance} a single number is returned.
#' 
#' \code{distance_table} returns a named list with two entries: \code{res} is
#' a numeric vector, the histogram of distances, \code{unconnected} is a
#' numeric scalar, the number of pairs for which the first vertex is not
#' reachable from the second. The sum of the two entries is always \eqn{n(n-1)}
#' for directed graphs and \eqn{n(n-1)/2} for undirected graphs.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @references West, D.B. (1996). \emph{Introduction to Graph Theory.} Upper
#' Saddle River, N.J.: Prentice Hall.
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- make_ring(10)
#' distances(g)
#' shortest_paths(g, 5)
#' all_shortest_paths(g, 1, 6:8)
#' mean_distance(g)
#' ## Weighted shortest paths
#' el <- matrix(nc=3, byrow=TRUE,
#'              c(1,2,0, 1,3,2, 1,4,1, 2,3,0, 2,5,5, 2,6,2, 3,2,1, 3,4,1,
#'                3,7,1, 4,3,0, 4,7,2, 5,6,2, 5,8,8, 6,3,2, 6,7,1, 6,9,1,
#'                6,10,3, 8,6,1, 8,9,1, 9,10,4) )
#' g2 <- add_edges(make_empty_graph(10), t(el[,1:2]), weight=el[,3])
#' distances(g2, mode="out")
#' 
distances <- function(graph, v=V(graph), to=V(graph),
                           mode=c("all", "out", "in"),
                           weights=NULL,
                           algorithm=c("automatic", "unweighted", "dijkstra",
                             "bellman-ford", "johnson")) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  v <- as.igraph.vs(graph, v)
  to <- as.igraph.vs(graph, to)
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)  
  algorithm <- igraph.match.arg(algorithm)
  algorithm <- switch(algorithm, "automatic"=0, "unweighted"=1,
                      "dijkstra"=2, "bellman-ford"=3, "johnson"=4)
  
  if (is.null(weights)) {
    if ("weight" %in% edge_attr_names(graph)) {
      weights <- as.numeric(E(graph)$weight)
    }
  } else {
    if (length(weights)==1 && is.na(weights)) {
      weights <- NULL
    } else {
      weights <- as.numeric(weights)
    }
  }

  if (! is.null(weights) && algorithm==1) {
    weights <- NULL
    warning("Unweighted algorithm chosen, weights ignored")
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_shortest_paths", graph, v-1, to-1,
               as.numeric(mode), weights, as.numeric(algorithm),
               PACKAGE="igraph")
  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    rownames(res) <- V(graph)$name[v]
    colnames(res) <- V(graph)$name[to]
  }
  res
}

#' @rdname distances
#' @param from Numeric constant, the vertex from or to the shortest paths will
#' be calculated. Note that right now this is not a vector of vertex ids, but
#' only a single vertex.
#' @param output Character scalar, defines how to report the shortest paths.
#' \dQuote{vpath} means that the vertices along the paths are reported, this
#' form was used prior to igraph version 0.6. \dQuote{epath} means that the
#' edges along the paths are reported. \dQuote{both} means that both forms are
#' returned, in a named list with components \dQuote{vpath} and \dQuote{epath}.
#' @param predecessors Logical scalar, whether to return the predecessor vertex
#' for each vertex. The predecessor of vertex \code{i} in the tree is the
#' vertex from which vertex \code{i} was reached. The predecessor of the start
#' vertex (in the \code{from} argument) is itself by definition. If the
#' predecessor is zero, it means that the given vertex was not reached from the
#' source during the search. Note that the search terminates if all the
#' vertices in \code{to} are reached.
#' @param inbound.edges Logical scalar, whether to return the inbound edge for
#' each vertex. The inbound edge of vertex \code{i} in the tree is the edge via
#' which vertex \code{i} was reached. The start vertex and vertices that were
#' not reached during the search will have zero in the corresponding entry of
#' the vector. Note that the search terminates if all the vertices in \code{to}
#' are reached.
#' @export

shortest_paths <- function(graph, from, to=V(graph),
                               mode=c("out", "all", "in"),
                               weights=NULL,
                               output=c("vpath", "epath", "both"),
                               predecessors=FALSE, inbound.edges=FALSE) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)
  output <- igraph.match.arg(output)
  output <- switch(output, "vpath"=0, "epath"=1, "both"=2)

  if (is.null(weights)) {
    if ("weight" %in% edge_attr_names(graph)) {
      weights <- as.numeric(E(graph)$weight)
    }
  } else {
    if (length(weights)==1 && is.na(weights)) {
      weights <- NULL
    } else {
      weights <- as.numeric(weights)
    }
  }
  
  to <- as.igraph.vs(graph, to)-1
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_get_shortest_paths", graph,
               as.igraph.vs(graph, from)-1, to, as.numeric(mode),
               as.numeric(length(to)), weights, as.numeric(output),
               as.logical(predecessors), as.logical(inbound.edges), 
               PACKAGE="igraph")

  if (!is.null(res$vpath)) {
    res$vpath <- lapply(res$vpath, function(x) x+1)
  }
  if (!is.null(res$epath)) {
    res$epath <- lapply(res$epath, function(x) x+1)
  }
  if (!is.null(res$predecessors)) {
    res$predecessors <- res$predecessors + 1
  }
  if (!is.null(res$inbound_edges)) {
    res$inbound_edges <- res$inbound_edges + 1
  }

  res
}

#' @export

all_shortest_paths <- function(graph, from,
                                   to=V(graph),
                                   mode=c("out", "all", "in"),
				   weights=NULL) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)

  if (is.null(weights)) {
    if ("weight" %in% edge_attr_names(graph)) {
      weights <- as.numeric(E(graph)$weight)
    }
  } else {
    if (length(weights)==1 && is.na(weights)) {
      weights <- NULL
    } else {
      weights <- as.numeric(weights)
    }
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  if (is.null(weights)) {
    res <- .Call("R_igraph_get_all_shortest_paths", graph,
                 as.igraph.vs(graph, from)-1, as.igraph.vs(graph, to)-1,
                 as.numeric(mode), PACKAGE="igraph")
  } else {
    res <- .Call("R_igraph_get_all_shortest_paths_dijkstra", graph, 
                 as.igraph.vs(graph, from)-1, as.igraph.vs(graph, to)-1,
                 weights, as.numeric(mode), PACKAGE="igraph")
  }       
  res
}

#' In- or out- component of a vertex
#' 
#' Finds all vertices reachable from a given vertex, or the opposite: all
#' vertices from which a given vertex is reachable via a directed path.
#' 
#' A breadh-first search is conducted starting from vertex \code{v}.
#' 
#' @aliases subcomponent
#' @param graph The graph to analyze.
#' @param v The vertex to start the search from.
#' @param mode Character string, either \dQuote{in}, \dQuote{out} or
#' \dQuote{all}. If \dQuote{in} all vertices from which \code{v} is reachable
#' are listed. If \dQuote{out} all vertices reachable from \code{v} are
#' returned. If \dQuote{all} returns the union of these. It is ignored for
#' undirected graphs.
#' @return Numeric vector, the ids of the vertices in the same component as
#' \code{v}.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{components}}
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- sample_gnp(100, 1/200)
#' subcomponent(g, 1, "in")
#' subcomponent(g, 1, "out")
#' subcomponent(g, 1, "all")

subcomponent <- function(graph, v, mode=c("all", "out", "in")) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_subcomponent", graph, as.igraph.vs(graph, v)-1,
               as.numeric(mode),
               PACKAGE="igraph")
  res+1
}



#' Subgraph of a graph
#' 
#' \code{subgraph} creates a subgraph of a graph, containing only the specified
#' vertices and all the edges among them.
#' 
#' \code{induced_subgraph} calculates the induced subgraph of a set of vertices
#' in a graph. This means that exactly the specified vertices and all the edges
#' between then will be kept in the result graph.
#' 
#' \code{subgraph.edges} calculates the subgraph of a graph. For this function
#' one can specify the vertices and edges to keep. This function will be
#' renamed to \code{subgraph} in the next major version of igraph.
#' 
#' The \code{subgraph} function does the same as \code{induced.graph} currently
#' (assuming \sQuote{\code{auto}} as the \code{impl} argument), but it is
#' deprecated and will be removed in the next major version of igraph.
#' 
#' @aliases subgraph induced.subgraph subgraph.edges induced_subgraph
#' @param graph The original graph.
#' @param v Numeric vector, the vertices of the original graph which will
#' form the subgraph.
#' @return A new graph object.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- make_ring(10)
#' g2 <- induced_subgraph(g, 1:7)
#' g3 <- subgraph.edges(g, 1:5, 1:5)
#' 
subgraph <- function(graph, v) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_subgraph", graph, as.igraph.vs(graph, v)-1,
        PACKAGE="igraph")
}

#' @rdname subgraph
#' @param vids Numeric vector, the vertices of the original graph which will
#' form the subgraph.
#' @param impl Character scalar, to choose between two implementation of the
#' subgraph calculation. \sQuote{\code{copy_and_delete}} copies the graph
#' first, and then deletes the vertices and edges that are not included in the
#' result graph. \sQuote{\code{create_from_scratch}} searches for all vertices
#' and edges that must be kept and then uses them to create the graph from
#' scratch. \sQuote{\code{auto}} chooses between the two implementations
#' automatically, using heuristics based on the size of the original and the
#' result graph.
#' @export

induced_subgraph <- function(graph, vids, impl=c("auto", "copy_and_delete", "create_from_scratch")) {
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  vids <- as.igraph.vs(graph, vids)
  impl <- switch(igraph.match.arg(impl), "auto"=0, "copy_and_delete"=1, "create_from_scratch"=2)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_induced_subgraph", graph, vids-1, impl,
        PACKAGE="igraph")

  res
}

#' @rdname subgraph
#' @param eids The edge ids of the edges that will be kept in the result graph.
#' @param delete.vertices Logical scalar, whether to remove vertices that do
#' not have any adjacent edges in \code{eids}.
#' @export

subgraph.edges <- function(graph, eids, delete.vertices=TRUE) {
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  eids <- as.igraph.es(graph, eids)
  delete.vertices <- as.logical(delete.vertices)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_subgraph_edges", graph, eids-1, delete.vertices,
        PACKAGE="igraph")

  res
}

#' @rdname betweenness
#' @param vids The vertices for which the vertex betweenness estimation will be
#' calculated.
#' @param cutoff The maximum path length to consider when calculating the
#' betweenness. If zero or negative then there is no such limit.
#' @export

estimate_betweenness <- function(graph, vids=V(graph), directed=TRUE, cutoff, weights=NULL, nobigint=TRUE) {
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  vids <- as.igraph.vs(graph, vids)
  directed <- as.logical(directed)
  cutoff <- as.numeric(cutoff)
  if (is.null(weights) && "weight" %in% edge_attr_names(graph)) { 
  weights <- E(graph)$weight 
  } 
  if (!is.null(weights) && any(!is.na(weights))) { 
  weights <- as.numeric(weights) 
  } else { 
  weights <- NULL 
  }
  nobigint <- as.logical(nobigint)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_betweenness_estimate", graph, vids-1, directed, cutoff, weights, nobigint,
        PACKAGE="igraph")
  if (igraph_opt("add.vertex.names") && is_named(graph)) { 
  names(res) <- vertex_attr(graph, "name", vids) 
  }
  res
}



#' Vertex and edge betweenness centrality
#' 
#' The vertex and edge betweenness are (roughly) defined by the number of
#' geodesics (shortest paths) going through a vertex or an edge.
#' 
#' The vertex betweenness of vertex \eqn{v}{\code{v}} is defined by
#' 
#' \deqn{\sum_{i\ne j, i\ne v, j\ne v} g_{ivj}/g_{ij}}{sum( g_ivj / g_ij,
#' i!=j,i!=v,j!=v)}
#' 
#' The edge betweenness of edge \eqn{e}{\code{e}} is defined by
#' 
#' \deqn{\sum_{i\ne j} g{iej}/g_{ij}.}{sum( g_iej / g_ij, i!=j).}
#' 
#' \code{betweenness} calculates vertex betweenness, \code{edge_betweenness}
#' calculates edge betweenness.
#' 
#' \code{estimate_betweenness} only considers paths of length \code{cutoff} or
#' smaller, this can be run for larger graphs, as the running time is not
#' quadratic (if \code{cutoff} is small). If \code{cutoff} is zero or negative
#' then the function calculates the exact betweenness scores.
#' 
#' \code{estimate_edge_betweenness} is similar, but for edges.
#' 
#' For calculating the betweenness a similar algorithm to the one proposed by
#' Brandes (see References) is used.
#' 
#' @aliases betweenness edge.betweenness betweenness.estimate
#' edge.betweenness.estimate edge_betweenness estimate_betweenness
#' estimate_edge_betweenness
#' @param graph The graph to analyze.
#' @param v The vertices for which the vertex betweenness will be calculated.
#' @param directed Logical, whether directed paths should be considered while
#' determining the shortest paths.
#' @param weights Optional positive weight vector for calculating weighted
#' betweenness. If the graph has a \code{weight} edge attribute, then this is
#' used by default.
#' @param nobigint Logical scalar, whether to use big integers during the
#' calculation. This is only required for lattice-like graphs that have very
#' many shortest paths between a pair of vertices. If \code{TRUE} (the
#' default), then big integers are not used.
#' @param normalized Logical scalar, whether to normalize the betweenness
#' scores. If \code{TRUE}, then the results are normalized according to
#' \deqn{B^n=\frac{2B}{n^2-3n+2}}{Bnorm=2*B/(n*n-3*n+2)}, where
#' \eqn{B^n}{Bnorm} is the normalized, \eqn{B} the raw betweenness, and \eqn{n}
#' is the number of vertices in the graph.
#' @return A numeric vector with the betweenness score for each vertex in
#' \code{v} for \code{betweenness}.
#' 
#' A numeric vector with the edge betweenness score for each edge in \code{e}
#' for \code{edge_betweenness}.
#' 
#' \code{estimate_betweenness} returns the estimated betweenness scores for
#' vertices in \code{vids}, \code{estimate_edge_betweenness} the estimated edge
#' betweenness score for \emph{all} edges; both in a numeric vector.
#' @note \code{edge_betweenness} might give false values for graphs with
#' multiple edges.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{closeness}}, \code{\link{degree}}
#' @references Freeman, L.C. (1979). Centrality in Social Networks I:
#' Conceptual Clarification. \emph{Social Networks}, 1, 215-239.
#' 
#' Ulrik Brandes, A Faster Algorithm for Betweenness Centrality. \emph{Journal
#' of Mathematical Sociology} 25(2):163-177, 2001.
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- sample_gnp(10, 3/10)
#' betweenness(g)
#' edge_betweenness(g)
#' 
betweenness <- function(graph, v=V(graph), directed=TRUE, weights=NULL,
                        nobigint=TRUE, normalized=FALSE) {
  
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  v <- as.igraph.vs(graph, v)
  if (is.null(weights) && "weight" %in% edge_attr_names(graph)) {
    weights <- E(graph)$weight
  }
  if (!is.null(weights) && any(!is.na(weights))) {
    weights <- as.numeric(weights)
  } else {
    weights <- NULL
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_betweenness", graph, v-1,
               as.logical(directed), weights, as.logical(nobigint),
               PACKAGE="igraph")
  if (normalized) {
    vc <- vcount(graph)
    if (is_directed(graph) && directed) {
      res <- res / ( vc*vc-3*vc+2)
    } else {
      res <- 2*res / ( vc*vc-3*vc+2)
    }
  }
  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    names(res) <- V(graph)$name[v]
  }
  res
}



#' Transitivity of a graph
#' 
#' Transitivity measures the probability that the adjacent vertices of a vertex
#' are connected. This is sometimes also called the clustering coefficient.
#' 
#' Note that there are essentially two classes of transitivity measures, one is
#' a vertex-level, the other a graph level property.
#' 
#' There are several generalizations of transitivity to weighted graphs, here
#' we use the definition by A. Barrat, this is a local vertex-level quantity,
#' its formula is
#' 
#' \deqn{C_i^w=\frac{1}{s_i(k_i-1)}\sum_{j,h}\frac{w_{ij}+w_{ih}}{2}a_{ij}a_{ih}a_{jh}}{
#' weighted C_i = 1/s_i 1/(k_i-1) sum( (w_ij+w_ih)/2 a_ij a_ih a_jh, j, h)}
#' 
#' \eqn{s_i}{s_i} is the strength of vertex \eqn{i}{i}, see
#' \code{\link{strength}}, \eqn{a_{ij}}{a_ij} are elements of the
#' adjacency matrix, \eqn{k_i}{k_i} is the vertex degree, \eqn{w_{ij}}{w_ij}
#' are the weights.
#' 
#' This formula gives back the normal not-weighted local transitivity if all
#' the edge weights are the same.
#' 
#' The \code{barrat} type of transitivity does not work for graphs with
#' multiple and/or loop edges. If you want to calculate it for a directed
#' graph, call \code{\link{as.undirected}} with the \code{collapse} mode first.
#' 
#' @param graph The graph to analyze.
#' @param type The type of the transitivity to calculate. Possible values:
#' \describe{ \item{list("global")}{The global transitivity of an undirected
#' graph (directed graphs are considered as undirected ones as well). This is
#' simply the ratio of the triangles and the connected triples in the graph.
#' For directed graph the direction of the edges is ignored. }
#' \item{list("local")}{The local transitivity of an undirected graph, this is
#' calculated for each vertex given in the \code{vids} argument. The local
#' transitivity of a vertex is the ratio of the triangles connected to the
#' vertex and the triples centered on the vertex. For directed graph the
#' direction of the edges is ignored. } \item{list("undirected")}{This is the
#' same as \code{global}.} \item{list("globalundirected")}{This is the same as
#' \code{global}.} \item{list("localundirected")}{This is the same as
#' \code{local}.} \item{list("barrat")}{The weighted transitivity as defined A.
#' Barrat. See details below.} \item{list("weighted")}{The same as
#' \code{barrat}.} }
#' @param vids The vertex ids for the local transitivity will be calculated.
#' This will be ignored for global transitivity types.  The default value is
#' \code{NULL}, in this case all vertices are considered. It is slightly faster
#' to supply \code{NULL} here than \code{V(graph)}.
#' @param weights Optional weights for weighted transitivity. It is ignored for
#' other transitivity measures. If it is \code{NULL} (the default) and the
#' graph has a \code{weight} edge attribute, then it is used automatically.
#' @param isolates Character scalar, defines how to treat vertices with degree
#' zero and one. If it is \sQuote{\code{NaN}} then they local transitivity is
#' reported as \code{NaN} and they are not included in the averaging, for the
#' transitivity types that calculate an average. If there are no vertices with
#' degree two or higher, then the averaging will still result \code{NaN}. If it
#' is \sQuote{\code{zero}}, then we report 0 transitivity for them, and they
#' are included in the averaging, if an average is calculated.
#' @return For \sQuote{\code{global}} a single number, or \code{NaN} if there
#' are no connected triples in the graph.
#' 
#' For \sQuote{\code{local}} a vector of transitivity scores, one for each
#' vertex in \sQuote{\code{vids}}.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @references Wasserman, S., and Faust, K. (1994). \emph{Social Network
#' Analysis: Methods and Applications.} Cambridge: Cambridge University Press.
#' 
#' Alain Barrat, Marc Barthelemy, Romualdo Pastor-Satorras, Alessandro
#' Vespignani: The architecture of complex weighted networks, Proc. Natl. Acad.
#' Sci. USA 101, 3747 (2004)
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- make_ring(10)
#' transitivity(g)
#' g2 <- sample_gnp(1000, 10/1000)
#' transitivity(g2)   # this is about 10/1000
#' 
#' # Weighted version, the figure from the Barrat paper
#' gw <- graph_from_literal(A-B:C:D:E, B-C:D, C-D)
#' E(gw)$weight <- 1
#' E(gw)[ V(gw)[name == "A"] %--% V(gw)[name == "E" ] ]$weight <- 5
#' transitivity(gw, vids="A", type="local")
#' transitivity(gw, vids="A", type="weighted")
#' 
#' # Weighted reduces to "local" if weights are the same
#' gw2 <- sample_gnp(1000, 10/1000)
#' E(gw2)$weight <- 1
#' t1 <- transitivity(gw2, type="local")
#' t2 <- transitivity(gw2, type="weighted")
#' all(is.na(t1) == is.na(t2))
#' all(na.omit(t1 == t2))
#' 
transitivity <- function(graph, type=c("undirected", "global", "globalundirected",
                                  "localundirected", "local", "average",
                                  "localaverage", "localaverageundirected",
                                  "barrat", "weighted"),
                         vids=NULL, weights=NULL, isolates=c("NaN", "zero")) {
  
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  type <- igraph.match.arg(type)
  type <- switch(type, "undirected"=0, "global"=0, "globalundirected"=0,
                 "localundirected"=1, "local"=1, "average"=2,
                 "localaverage"=2, "localaverageundirected"=2, "barrat"=3,
                 "weighted"=3)

  if (is.null(weights) && "weight" %in% edge_attr_names(graph)) {
    weights <- E(graph)$weight
  }
  if (!is.null(weights) && any(!is.na(weights))) {
    weights <- as.numeric(weights)
  } else {
    weights <- NULL
  }

  isolates <- igraph.match.arg(isolates)
  isolates <- as.double(switch(isolates, "nan"=0, "zero"=1))

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  if (type==0) {
    .Call("R_igraph_transitivity_undirected", graph, isolates,
          PACKAGE="igraph")
  } else if (type==1) {
    if (is.null(vids)) {
      .Call("R_igraph_transitivity_local_undirected_all", graph, isolates,
            PACKAGE="igraph")
    } else {
      vids <- as.igraph.vs(graph, vids)-1
      .Call("R_igraph_transitivity_local_undirected", graph, vids,
            isolates, PACKAGE="igraph")
    }
  } else if (type==2) {
    .Call("R_igraph_transitivity_avglocal_undirected", graph, isolates,
          PACKAGE="igraph")
  } else if (type==3) {
    if (is.null(vids)) { vids <- V(graph) }
    vids <- as.igraph.vs(graph, vids)-1
    if (is.null(weights)) {
      .Call("R_igraph_transitivity_local_undirected", graph, vids,
            isolates, PACKAGE="igraph")
    } else { 
      .Call("R_igraph_transitivity_barrat", graph, vids, weights,
            isolates, PACKAGE="igraph")
    }
  }
}

## Generated by stimulus now
## laplacian_matrix <- function(graph, normalized=FALSE) {

##   if (!is_igraph(graph)) {
##     stop("Not a graph object")
##   }
  
##   on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
##   .Call("R_igraph_laplacian", graph, as.logical(normalized),
##         PACKAGE="igraph")
## }
  
## OLD implementation
## laplacian_matrix <- function(graph, normalized=FALSE) {

##   if (!is_igraph(graph)) {
##     stop("Not a graph object")
##   }
##   if (is_directed(graph)) {
##     warning("Laplacian of a directed graph???")
##   }

##   M <- as_adj(graph)
##   if (!normalized) {
##     M <- structure(ifelse(M>0, -1, 0), dim=dim(M))
##     diag(M) <- degree(graph)
##   } else {
##     deg <- degree(graph)
##     deg <- outer(deg, deg, "*")
##     M <- structure(ifelse(M>0, -1/deg, 0))
##     diag(M) <- 1
##   }
  
##   M
## }

## Structural holes a'la Burt, code contributed by
## Jeroen Bruggeman

## constraint.orig <- function(graph, nodes=V(graph), attr=NULL) {

##   if (!is_igraph(graph)) {
##     stop("Not a graph object")
##   }

##   idx <- degree(graph) != 0
##   A <- as_adj(graph, attr=attr)
##   A <- A[idx, idx]
##   n <- sum(idx)
  
##   one <- c(rep(1,n))
##   CZ <- A + t(A)
##   cs <- CZ %*% one                      # degree of vertices
##   ics <- 1/cs
##   CS <- ics %*% t(one)                  # 1/degree of vertices
##   P <- CZ * CS  #intermediate result: proportionate tie strengths
##   PSQ <- P%*%P #sum paths of length two
##   P.bi <- as.numeric(P>0)  #exclude paths to non-contacts (& reflexive):
##   PC <- (P + (PSQ*P.bi))^2  #dyadic constraint
##   ci <- PC %*% one      #overall constraint
##   dim(ci) <- NULL

##   ci2 <- numeric(vcount(graph))
##   ci2[idx] <- ci
##   ci2[!idx] <- NaN
##   ci2[nodes+1]
## }

## Newest implementation, hopefully correct, there is a C implementation
## now so we don't need this

## constraint.old <- function(graph, nodes=V(graph)) {

##   if (!is_igraph(graph)) {
##     stop("Not a graph object")
##   }

##   nodes <- as.numeric(nodes)
##   res <- numeric(length(nodes))
##   deg <- degree(graph, mode="all", loops=FALSE)

##   not <- function(i, v) v[ v!=i ]

##   for (a in seq(along=nodes)) {
##     i <- nodes[a]
    
##     first <- not(i, neighbors(graph, i, mode="all"))
##     first <- unique(first)
##     for (b in seq(along=first)) {
##       j <- first[b]

##       ## cj is the contribution of j
##       cj <- is_connected_to(graph, i, j)      / deg[i+1]
##       cj <- cj + is_connected_to(graph, j, i) / deg[i+1]

##       second <- not(i, not(j, neighbors(graph, j, mode="all")))
##       for (c in seq(along=second)) {
##         q <- second[c]
##         cj <- cj + is_connected_to(graph, i, q) / deg[q+1] / deg[i+1]
##         cj <- cj + is_connected_to(graph, q, i) / deg[q+1] / deg[i+1]
##       }
                            
##       ## Ok, we have the total contribution of j
##       res[a] <- res[a] + cj*cj
##     }
##   }

##   if (!is_directed(graph)) {
##     res <- res/4
##   }
##   res
## }



#' Burt's constraint
#' 
#' Given a graph, \code{constraint} calculates Burt's constraint for each
#' vertex.
#'
#' Burt's constraint is higher if ego has less, or mutually
#' stronger related (i.e. more redundant) contacts. Burt's measure of
#' constraint, \eqn{C_i}{C[i]}, of vertex \eqn{i}'s ego network
#' \eqn{V_i}{V[i]}, is defined for directed and valued graphs,
#' \deqn{C_i=\sum_{j \in V_i \setminus \{i\}} (p_{ij}+\sum_{q \in V_i
#'     \setminus \{i,j\}} p_{iq} p_{qj})^2}{
#'   C[i] = sum( [sum( p[i,j] + p[i,q] p[q,j], q in V[i], q != i,j )]^2, j in
#'   V[i], j != i).
#' }
#' for a graph of order (ie. number of vertices) \eqn{N}, where
#' proportional tie strengths are defined as 
#' \deqn{p_{ij} = \frac{a_{ij}+a_{ji}}{\sum_{k \in V_i \setminus \{i\}}(a_{ik}+a_{ki})},}{
#'   p[i,j]=(a[i,j]+a[j,i]) / sum(a[i,k]+a[k,i], k in V[i], k != i),
#' }
#' \eqn{a_{ij}}{a[i,j]} are elements of \eqn{A} and the latter being the
#' graph adjacency matrix. For isolated vertices, constraint is undefined.
#' 
#' @param graph A graph object, the input graph.
#' @param nodes The vertices for which the constraint will be calculated.
#' Defaults to all vertices.
#' @param weights The weights of the edges. If this is \code{NULL} and there is
#' a \code{weight} edge attribute this is used. If there is no such edge
#' attribute all edges will have the same weight.
#' @return A numeric vector of constraint scores
#' @author Jeroen Bruggeman
#' (\url{https://sites.google.com/site/jebrug/jeroen-bruggeman-social-science})
#' and Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @references Burt, R.S. (2004). Structural holes and good ideas.
#' \emph{American Journal of Sociology} 110, 349-399.
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- sample_gnp(20, 5/20)
#' constraint(g)
#' 
constraint <- function(graph, nodes=V(graph), weights=NULL) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  nodes <- as.igraph.vs(graph, nodes)
  
  if (is.null(weights)) {
    if ("weight" %in% edge_attr_names(graph)) {
      weights <- E(graph)$weight
    }
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_constraint", graph, nodes-1, as.numeric(weights),
               PACKAGE="igraph")
  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    names(res) <- V(graph)$name[nodes]
  }
  res
}



#' Reciprocity of graphs
#' 
#' Calculates the reciprocity of a directed graph.
#' 
#' The measure of reciprocity defines the proporsion of mutual connections, in
#' a directed graph. It is most commonly defined as the probability that the
#' opposite counterpart of a directed edge is also included in the graph. Or in
#' adjacency matrix notation: \eqn{\sum_{ij} (A\cdot A')_{ij}}{sum(i, j,
#' (A.*A')ij) / sum(i, j, Aij)}, where \eqn{A\cdot A'}{A.*A'} is the
#' element-wise product of matrix \eqn{A} and its transpose. This measure is
#' calculated if the \code{mode} argument is \code{default}.
#' 
#' Prior to igraph version 0.6, another measure was implemented, defined as the
#' probability of mutual connection between a vertex pair, if we know that
#' there is a (possibly non-mutual) connection between them. In other words,
#' (unordered) vertex pairs are classified into three groups: (1)
#' not-connected, (2) non-reciprocaly connected, (3) reciprocally connected.
#' The result is the size of group (3), divided by the sum of group sizes
#' (2)+(3). This measure is calculated if \code{mode} is \code{ratio}.
#' 
#' @param graph The graph object.
#' @param ignore.loops Logical constant, whether to ignore loop edges.
#' @param mode See below.
#' @return A numeric scalar between zero and one.
#' @author Tamas Nepusz \email{ntamas@@gmail.com} and Gabor Csardi
#' \email{csardi.gabor@@gmail.com}
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- sample_gnp(20, 5/20, directed=TRUE)
#' reciprocity(g)
#' 
reciprocity <- function(graph, ignore.loops=TRUE,
                        mode=c("default", "ratio")) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- switch(igraph.match.arg(mode), 'default'=0, 'ratio'=1)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_reciprocity", graph, as.logical(ignore.loops),
        as.numeric(mode), PACKAGE="igraph")
}


bonpow.dense <- function(graph, nodes=V(graph),
                         loops=FALSE, exponent=1,
                         rescale=FALSE, tol=1e-7){

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }  
  
  d <- as_adj(graph)
  if (!loops) {
    diag(d) <- 0
  }
  n <- vcount(graph)
  id <- matrix(0,nrow=n,ncol=n)
  diag(id) <- 1

#  ev <- apply(solve(id-exponent*d,tol=tol)%*%d,1,sum)
  ev <- solve(id-exponent*d, tol=tol) %*% apply(d,1,sum)
  if(rescale) {
    ev <- ev/sum(ev)
  } else {
    ev <- ev*sqrt(n/sum((ev)^2))
  } 
  ev[as.numeric(nodes)]
}

bonpow.sparse <- function(graph, nodes=V(graph), loops=FALSE,
                          exponent=1, rescale=FALSE, tol=1e-07) {

  ## remove loops if requested
  if (!loops) {
    graph <- simplify(graph, remove.multiple=FALSE, remove.loops=TRUE)
  }

  vg <- vcount(graph)
  
  ## sparse adjacency matrix
  d <- as_adj(graph, sparse=TRUE)

  ## sparse identity matrix
  id <- Matrix::Diagonal(vg)

  ## solve it
  ev <- Matrix::solve(id - exponent * d, degree(graph, mode="out"), tol=tol)

  if (rescale) {
    ev <- ev/sum(ev)
  } else {
    ev <- ev * sqrt(vcount(graph)/sum((ev)^2))
  }

  ev[as.numeric(nodes)]
}



#' Find Bonacich Power Centrality Scores of Network Positions
#' 
#' \code{power_centrality} takes a graph (\code{dat}) and returns the Boncich power
#' centralities of positions (selected by \code{nodes}).  The decay rate for
#' power contributions is specified by \code{exponent} (1 by default).
#' 
#' Bonacich's power centrality measure is defined by
#' \eqn{C_{BP}\left(\alpha,\beta\right)=\alpha\left(\mathbf{I}-\beta\mathbf{A}\right)^{-1}\mathbf{A}\mathbf{1}}{C_BP(alpha,beta)=alpha
#' (I-beta A)^-1 A 1}, where \eqn{\beta}{beta} is an attenuation parameter (set
#' here by \code{exponent}) and \eqn{\mathbf{A}}{A} is the graph adjacency
#' matrix.  (The coefficient \eqn{\alpha}{alpha} acts as a scaling parameter,
#' and is set here (following Bonacich (1987)) such that the sum of squared
#' scores is equal to the number of vertices.  This allows 1 to be used as a
#' reference value for the ``middle'' of the centrality range.)  When
#' \eqn{\beta \rightarrow }{beta->1/lambda_A1}\eqn{
#' 1/\lambda_{\mathbf{A}1}}{beta->1/lambda_A1} (the reciprocal of the largest
#' eigenvalue of \eqn{\mathbf{A}}{A}), this is to within a constant multiple of
#' the familiar eigenvector centrality score; for other values of \eqn{\beta},
#' the behavior of the measure is quite different.  In particular, \eqn{\beta}
#' gives positive and negative weight to even and odd walks, respectively, as
#' can be seen from the series expansion
#' \eqn{C_{BP}\left(\alpha,\beta\right)=\alpha \sum_{k=0}^\infty \beta^k
#' }{C_BP(alpha,beta) = alpha sum( beta^k A^(k+1) 1, k in 0..infinity )}\eqn{
#' \mathbf{A}^{k+1} \mathbf{1}}{C_BP(alpha,beta) = alpha sum( beta^k A^(k+1) 1,
#' k in 0..infinity )} which converges so long as \eqn{|\beta|
#' }{|beta|<1/lambda_A1}\eqn{ < 1/\lambda_{\mathbf{A}1}}{|beta|<1/lambda_A1}.
#' The magnitude of \eqn{\beta}{beta} controls the influence of distant actors
#' on ego's centrality score, with larger magnitudes indicating slower rates of
#' decay.  (High rates, hence, imply a greater sensitivity to edge effects.)
#' 
#' Interpretively, the Bonacich power measure corresponds to the notion that
#' the power of a vertex is recursively defined by the sum of the power of its
#' alters.  The nature of the recursion involved is then controlled by the
#' power exponent: positive values imply that vertices become more powerful as
#' their alters become more powerful (as occurs in cooperative relations),
#' while negative values imply that vertices become more powerful only as their
#' alters become \emph{weaker} (as occurs in competitive or antagonistic
#' relations).  The magnitude of the exponent indicates the tendency of the
#' effect to decay across long walks; higher magnitudes imply slower decay.
#' One interesting feature of this measure is its relative instability to
#' changes in exponent magnitude (particularly in the negative case).  If your
#' theory motivates use of this measure, you should be very careful to choose a
#' decay parameter on a non-ad hoc basis.
#'
#' @aliases bonpow
#' @param graph the input graph.
#' @param nodes vertex sequence indicating which vertices are to be included in
#' the calculation.  By default, all vertices are included.
#' @param loops boolean indicating whether or not the diagonal should be
#' treated as valid data.  Set this true if and only if the data can contain
#' loops.  \code{loops} is \code{FALSE} by default.
#' @param exponent exponent (decay rate) for the Bonacich power centrality
#' score; can be negative
#' @param rescale if true, centrality scores are rescaled such that they sum to
#' 1.
#' @param tol tolerance for near-singularities during matrix inversion (see
#' \code{\link{solve}})
#' @param sparse Logical scalar, whether to use sparse matrices for the
#' calculation. The \sQuote{Matrix} package is required for sparse matrix
#' support
#' @return A vector, containing the centrality scores.
#' @note This function was ported (ie. copied) from the SNA package.
#' @section Warning : Singular adjacency matrices cause no end of headaches for
#' this algorithm; thus, the routine may fail in certain cases.  This will be
#' fixed when I get a better algorithm.  \code{power_centrality} will not symmetrize your
#' data before extracting eigenvectors; don't send this routine asymmetric
#' matrices unless you really mean to do so.
#' @author Carter T. Butts
#' (\url{http://www.faculty.uci.edu/profile.cfm?faculty_id=5057}), ported to
#' igraph by Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{eigen_centrality}} and \code{\link{alpha_centrality}}
#' @references Bonacich, P.  (1972).  ``Factoring and Weighting Approaches to
#' Status Scores and Clique Identification.'' \emph{Journal of Mathematical
#' Sociology}, 2, 113-120.
#' 
#' Bonacich, P.  (1987).  ``Power and Centrality: A Family of Measures.''
#' \emph{American Journal of Sociology}, 92, 1170-1182.
#' @keywords graphs
#' @export
#' @examples
#' 
#' # Generate some test data from Bonacich, 1987:
#' g.c <- graph( c(1,2,1,3,2,4,3,5), dir=FALSE)
#' g.d <- graph( c(1,2,1,3,1,4,2,5,3,6,4,7), dir=FALSE)
#' g.e <- graph( c(1,2,1,3,1,4,2,5,2,6,3,7,3,8,4,9,4,10), dir=FALSE)
#' g.f <- graph( c(1,2,1,3,1,4,2,5,2,6,2,7,3,8,3,9,3,10,4,11,4,12,4,13), dir=FALSE)
#' # Compute power centrality scores
#' for (e in seq(-0.5,.5, by=0.1)) {
#'   print(round(power_centrality(g.c, exp=e)[c(1,2,4)], 2))
#' }
#' 
#' for (e in seq(-0.4,.4, by=0.1)) {
#'   print(round(power_centrality(g.d, exp=e)[c(1,2,5)], 2))
#' }
#' 
#' for (e in seq(-0.4,.4, by=0.1)) {
#'   print(round(power_centrality(g.e, exp=e)[c(1,2,5)], 2))
#' }
#' 
#' for (e in seq(-0.4,.4, by=0.1)) {
#'   print(round(power_centrality(g.f, exp=e)[c(1,2,5)], 2))
#' }
#' 
power_centrality <- function(graph, nodes=V(graph),
                   loops=FALSE, exponent=1,
                   rescale=FALSE, tol=1e-7, sparse=TRUE){

  nodes <- as.igraph.vs(graph, nodes)
  if (sparse) {
    res <- bonpow.sparse(graph, nodes, loops, exponent, rescale, tol)
  }  else {
    res <- bonpow.dense(graph, nodes, loops, exponent, rescale, tol)
  }

  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    names(res) <- vertex_attr(graph, "name", nodes)
  }
  
  res
}

alpha.centrality.dense <- function(graph, nodes=V(graph), alpha=1,
                                   loops=FALSE, exo=1, weights=NULL,
                                   tol=1e-7) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  exo <- rep(exo, length=vcount(graph))
  exo <- matrix(exo, ncol=1)

  if (is.null(weights) && "weight" %in% edge_attr_names(graph)) {
    ## weights == NULL and there is a "weight" edge attribute
    attr <- "weight"
  } else if (is.null(weights)) {
    ## weights == NULL, but there is no "weight" edge attribute
    attr <- NULL
  } else if (is.character(weights) && length(weights)==1) {
    ## name of an edge attribute, nothing to do
    attr <- "weight"
  } else if (any(!is.na(weights))) {
    ## weights != NULL and weights != rep(NA, x)
    graph <- set_edge_attr(graph, "weight", value=as.numeric(weights))
    attr <- "weight"
  } else {
    ## weights != NULL, but weights == rep(NA, x)
    attr <- NULL
  }

  d <- t(as_adj(graph, attr=attr, sparse=FALSE))
  if (!loops) {
    diag(d) <- 0
  }
  n <- vcount(graph)
  id <- matrix(0, nrow=n, ncol=n)
  diag(id) <- 1
  
  ev <- solve(id-alpha*d, tol=tol) %*% exo
  ev[as.numeric(nodes)]
}

alpha.centrality.sparse <- function(graph, nodes=V(graph), alpha=1,
                                   loops=FALSE, exo=1, weights=NULL,
                                   tol=1e-7) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  vc <- vcount(graph)

  if (!loops) {
    graph <- simplify(graph, remove.multiple=FALSE, remove.loops=TRUE)
  }  

  if (is.null(weights) && "weight" %in% edge_attr_names(graph)) {
    ## weights == NULL and there is a "weight" edge attribute
    attr <- "weight"
  } else if (is.null(weights)) {
    ## weights == NULL, but there is no "weight" edge attribute
    attr <- NULL
  } else if (is.character(weights) && length(weights)==1) {
    ## name of an edge attribute, nothing to do
    attr <- "weight"
  } else if (any(!is.na(weights))) {
    ## weights != NULL and weights != rep(NA, x)
    graph <- set_edge_attr(graph, "weight", value=as.numeric(weights))
    attr <- "weight"
  } else {
    ## weights != NULL, but weights == rep(NA, x)
    attr <- NULL
  }

  M <- Matrix::t(as_adj(graph, attr = attr, sparse = TRUE))
  M <- as(M, "dgCMatrix")
  
  ## Create an identity matrix
  M2 <- Matrix::sparseMatrix(dims=c(vc, vc), i=1:vc, j=1:vc, x=rep(1, vc))
  M2 <- as(M2, "dgCMatrix")

  ## exo
  exo <- cbind(rep(exo, length=vc))

  ## Solve the equation
  M3 <- M2-alpha*M
  r <- Matrix::solve(M3, tol=tol, exo)
  
  r[ as.numeric(nodes)]
}



#' Find Bonacich alpha centrality scores of network positions
#' 
#' \code{alpha_centrality} calculates the alpha centrality of some (or all)
#' vertices in a graph.
#' 
#' The alpha centrality measure can be considered as a generalization of
#' eigenvector centerality to directed graphs. It was proposed by Bonacich in
#' 2001 (see reference below).
#' 
#' The alpha centrality of the vertices in a graph is defined as the solution
#' of the following matrix equation: \deqn{x=\alpha A^T x+e,}{x=alpha t(A)x+e,}
#' where \eqn{A}{A} is the (not neccessarily symmetric) adjacency matrix of the
#' graph, \eqn{e}{e} is the vector of exogenous sources of status of the
#' vertices and \eqn{\alpha}{alpha} is the relative importance of the
#' endogenous versus exogenous factors.
#'
#' @aliases alpha.centrality
#' @param graph The input graph, can be directed or undirected
#' @param nodes Vertex sequence, the vertices for which the alpha centrality
#' values are returned. (For technical reasons they will be calculated for all
#' vertices, anyway.)
#' @param alpha Parameter specifying the relative importance of endogenous
#' versus exogenous factors in the determination of centrality. See details
#' below.
#' @param loops Whether to eliminate loop edges from the graph before the
#' calculation.
#' @param exo The exogenous factors, in most cases this is either a constant --
#' the same factor for every node, or a vector giving the factor for every
#' vertex. Note that too long vectors will be truncated and too short vectors
#' will be replicated to match the number of vertices.
#' @param weights A character scalar that gives the name of the edge attribute
#' to use in the adjacency matrix. If it is \code{NULL}, then the
#' \sQuote{weight} edge attribute of the graph is used, if there is one.
#' Otherwise, or if it is \code{NA}, then the calculation uses the standard
#' adjacency matrix.
#' @param tol Tolerance for near-singularities during matrix inversion, see
#' \code{\link{solve}}.
#' @param sparse Logical scalar, whether to use sparse matrices for the
#' calculation. The \sQuote{Matrix} package is required for sparse matrix
#' support
#' @return A numeric vector contaning the centrality scores for the selected
#' vertices.
#' @section Warning: Singular adjacency matrices cause problems for this
#' algorithm, the routine may fail is certain cases.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{eigen_centrality}} and \code{\link{power_centrality}}
#' @references Bonacich, P. and Paulette, L. (2001). ``Eigenvector-like
#' measures of centrality for asymmetric relations'' \emph{Social Networks},
#' 23, 191-201.
#' @export
#' @keywords graphs
#' @examples
#' 
#' # The examples from Bonacich's paper
#' g.1 <- graph( c(1,3,2,3,3,4,4,5) )
#' g.2 <- graph( c(2,1,3,1,4,1,5,1) )
#' g.3 <- graph( c(1,2,2,3,3,4,4,1,5,1) )
#' alpha_centrality(g.1)
#' alpha_centrality(g.2)
#' alpha_centrality(g.3,alpha=0.5)
#' 
alpha_centrality <- function(graph, nodes=V(graph), alpha=1,
                             loops=FALSE, exo=1, weights=NULL,
                             tol=1e-7, sparse=TRUE) {

  nodes <- as.igraph.vs(graph, nodes)
  if (sparse) {
    res <- alpha.centrality.sparse(graph, nodes, alpha, loops,
                                   exo, weights, tol)
  } else {
    res <- alpha.centrality.dense(graph, nodes, alpha, loops,
                                  exo, weights, tol)
  }
  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    names(res) <- vertex_attr(graph, "name", nodes)
  }
  res
}




#' Graph density
#' 
#' The density of a graph is the ratio of the number of edges and the number of
#' possible edges.
#' 
#' Note that this function may return strange results for graph with multiple
#' edges, density is ill-defined for graphs with multiple edges.
#'
#' @aliases graph.density
#' @param graph The input graph.
#' @param loops Logical constant, whether to allow loop edges in the graph. If
#' this is TRUE then self loops are considered to be possible. If this is FALSE
#' then we assume that the graph does not contain any loop edges and that loop
#' edges are not meaningful.
#' @return A real constant. This function returns \code{NaN} (=0.0/0.0) for an
#' empty graph with zero vertices.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{vcount}}, \code{\link{ecount}}, \code{\link{simplify}}
#' to get rid of the multiple and/or loop edges.
#' @references Wasserman, S., and Faust, K.  (1994).  Social Network Analysis:
#' Methods and Applications.  Cambridge: Cambridge University Press.
#' @export
#' @keywords graphs
#' @examples
#' 
#' g1 <- make_empty_graph(n=10)
#' g2 <- make_full_graph(n=10)
#' g3 <- sample_gnp(n=10, 0.4)
#' 
#' # loop edges
#' g <- graph( c(1,2, 2,2, 2,3) )
#' density(g, loops=FALSE)              # this is wrong!!!
#' density(g, loops=TRUE)               # this is right!!!
#' density(simplify(g), loops=FALSE)    # this is also right, but different
#' 
density <- function(graph, loops=FALSE) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }  
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_density", graph, as.logical(loops),
        PACKAGE="igraph")
}

#' @export

ego_size <- function(graph, order, nodes=V(graph),
                              mode=c("all", "out", "in"), mindist=0) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)
  mindist <- as.integer(mindist)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_neighborhood_size", graph, 
        as.igraph.vs(graph, nodes)-1, as.numeric(order), as.numeric(mode),
        mindist, PACKAGE="igraph")
}



#' Neighborhood of graph vertices
#' 
#' These functions find the vertices not farther than a given limit from
#' another fixed vertex, these are called the neighborhood of the vertex.
#' 
#' The neighborhood of a given order \code{o} of a vertex \code{v} includes all
#' vertices which are closer to \code{v} than the order. Ie. order 0 is always
#' \code{v} itself, order 1 is \code{v} plus its immediate neighbors, order 2
#' is order 1 plus the immediate neighbors of the vertices in order 1, etc.
#' 
#' \code{ego_size} calculates the size of the neighborhoods for the
#' given vertices with the given order.
#' 
#' \code{ego} calculates the neighborhoods of the given vertices with
#' the given order parameter.
#' 
#' \code{make_ego_graph} is creates (sub)graphs from all neighborhoods of
#' the given vertices with the given order parameter. This function preserves
#' the vertex, edge and graph attributes.
#' 
#' \code{connect} creates a new graph by connecting each vertex to
#' all other vertices in its neighborhood.
#' 
#' @aliases neighborhood neighborhood.size graph.neighborhood ego_graph
#' connect.neighborhood connect ego_size ego
#' @param graph The input graph.
#' @param order Integer giving the order of the neighborhood.
#' @param nodes The vertices for which the calculation is performed.
#' @param mode Character constatnt, it specifies how to use the direction of
#' the edges if a directed graph is analyzed. For \sQuote{out} only the
#' outgoing edges are followed, so all vertices reachable from the source
#' vertex in at most \code{order} steps are counted. For \sQuote{"in"} all
#' vertices from which the source vertex is reachable in at most \code{order}
#' steps are counted. \sQuote{"all"} ignores the direction of the edges. This
#' argument is ignored for undirected graphs.
#' @param mindist The minimum distance to include the vertex in the result.
#' @return \code{ego_size} returns with an integer vector.
#' 
#' \code{ego} returns with a list of integer vectors.
#' 
#' \code{make_ego_graph} returns with a list of graphs.
#' 
#' \code{connect} returns with a new graph object.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}, the first version was
#' done by Vincent Matossian
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- make_ring(10)
#' ego_size(g, 0, 1:3)
#' ego_size(g, 1, 1:3)
#' ego_size(g, 2, 1:3)
#' ego(g, 0, 1:3)
#' ego(g, 1, 1:3)
#' ego(g, 2, 1:3)
#' 
#' # attributes are preserved
#' V(g)$name <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j")
#' make_ego_graph(g, 2, 1:3)
#' 
#' # connecting to the neighborhood
#' g <- make_ring(10)
#' g <- connect(g, 2)
#' 
ego <- function(graph, order, nodes=V(graph),
                         mode=c("all", "out", "in"), mindist=0) {
  
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)
  mindist <- as.integer(mindist)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_neighborhood", graph, 
               as.igraph.vs(graph, nodes)-1, as.numeric(order),
               as.numeric(mode), mindist,
               PACKAGE="igraph")
  res <- lapply(res, function(x) x+1)
  res
}

#' @rdname ego
#' @export

make_ego_graph <- function(graph, order, nodes=V(graph),
                  mode=c("all", "out", "in"), mindist=0) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)
  mindist <- as.integer(mindist)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_neighborhood_graphs", graph, 
               as.igraph.vs(graph, nodes)-1, as.numeric(order),
               as.numeric(mode), mindist,
               PACKAGE="igraph")
  res
}



#' K-core decomposition of graphs
#' 
#' The k-core of graph is a maximal subgraph in which each vertex has at least
#' degree k. The coreness of a vertex is k if it belongs to the k-core but not
#' to the (k+1)-core.
#' 
#' The k-core of a graph is the maximal subgraph in which every vertex has at
#' least degree k. The cores of a graph form layers: the (k+1)-core is always a
#' subgraph of the k-core.
#' 
#' This function calculates the coreness for each vertex.
#'
#' @aliases graph.coreness
#' @param graph The input graph, it can be directed or undirected
#' @param mode The type of the core in directed graphs. Character constant,
#' possible values: \code{in}: in-cores are computed, \code{out}: out-cores are
#' computed, \code{all}: the corresponding undirected graph is considered. This
#' argument is ignored for undirected graphs.
#' @return Numeric vector of integer numbers giving the coreness of each
#' vertex.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{degree}}
#' @references Vladimir Batagelj, Matjaz Zaversnik: An O(m) Algorithm for Cores
#' Decomposition of Networks, 2002
#' 
#' Seidman S. B. (1983) Network structure and minimum degree, \emph{Social
#' Networks}, 5, 269--287.
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- make_ring(10)
#' g <- add_edges(g, c(1,2, 2,3, 1,3))
#' coreness(g) 		# small core triangle in a ring
#' 
coreness <- function(graph, mode=c("all", "out", "in")) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_coreness", graph, as.numeric(mode),
               PACKAGE="igraph")
  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    names(res) <- vertex_attr(graph, "name")
  }
  res
}



#' Topological sorting of vertices in a graph
#' 
#' A topological sorting of a directed acyclic graph is a linear ordering of
#' its nodes where each node comes before all nodes to which it has edges.
#' 
#' Every DAG has at least one topological sort, and may have many.  This
#' function returns a possible topological sort among them. If the graph is not
#' acyclic (it has at least one cycle), a partial topological sort is returned
#' and a warning is issued.
#'
#' @aliases topological.sort
#' @param graph The input graph, should be directed
#' @param mode Specifies how to use the direction of the edges.  For
#' \dQuote{\code{out}}, the sorting order ensures that each node comes before
#' all nodes to which it has edges, so nodes with no incoming edges go first.
#' For \dQuote{\code{in}}, it is quite the opposite: each node comes before all
#' nodes from which it receives edges. Nodes with no outgoing edges go first.
#' @return A numeric vector containing vertex ids in topologically sorted
#' order.
#' @author Tamas Nepusz \email{ntamas@@gmail.com} and Gabor Csardi
#' \email{csardi.gabor@@gmail.com} for the R interface
#' @keywords graphs
#' @export
#' @examples
#' 
#' g <- barabasi.game(100)
#' topo_sort(g)
#' 
topo_sort <- function(graph, mode=c("out", "all", "in")) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=1, "in"=2, "all"=3)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_topological_sorting", graph, as.numeric(mode),
               PACKAGE="igraph")
  res+1
}



#' Girth of a graph
#' 
#' The girth of a graph is the length of the shortest circle in it.
#' 
#' The current implementation works for undirected graphs only, directed graphs
#' are treated as undirected graphs. Loop edges and multiple edges are ignored.
#' If the graph is a forest (ie. acyclic), then zero is returned.
#' 
#' This implementation is based on Alon Itai and Michael Rodeh: Finding a
#' minimum circuit in a graph \emph{Proceedings of the ninth annual ACM
#' symposium on Theory of computing}, 1-10, 1977. The first implementation of
#' this function was done by Keith Briggs, thanks Keith.
#' 
#' @param graph The input graph. It may be directed, but the algorithm searches
#' for undirected circles anyway.
#' @param circle Logical scalar, whether to return the shortest circle itself.
#' @return A named list with two components: \item{girth}{Integer constant, the
#' girth of the graph, or 0 if the graph is acyclic.} \item{circle}{Numeric
#' vector with the vertex ids in the shortest circle.}
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @references Alon Itai and Michael Rodeh: Finding a minimum circuit in a
#' graph \emph{Proceedings of the ninth annual ACM symposium on Theory of
#' computing}, 1-10, 1977
#' @export
#' @keywords graphs
#' @examples
#' 
#' # No circle in a tree
#' g <- make_tree(1000, 3)
#' girth(g)
#' 
#' # The worst case running time is for a ring
#' g <- make_ring(100)
#' girth(g)
#' 
#' # What about a random graph?
#' g <- sample_gnp(1000, 1/1000)
#' girth(g)
#' 
girth <- function(graph, circle=TRUE) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_girth", graph, as.logical(circle),
        PACKAGE="igraph")
}

#' @export

which_loop <- function(graph, eids=E(graph)) {

  if (!is_igraph(graph)) {
    stop("Not a graph object");
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_is_loop", graph, as.igraph.es(graph, eids)-1,
        PACKAGE="igraph")
}



#' Find the multiple or loop edges in a graph
#' 
#' A loop edge is an edge from a vertex to itself. An edge is a multiple edge
#' if it has exactly the same head and tail vertices as another edge. A graph
#' without multiple and loop edges is called a simple graph.
#' 
#' \code{which_loop} decides whether the edges of the graph are loop edges.
#' 
#' \code{any_multiple} decides whether the graph has any multiple edges.
#' 
#' \code{which_multiple} decides whether the edges of the graph are multiple
#' edges.
#' 
#' \code{count_multiple} counts the multiplicity of each edge of a graph.
#' 
#' Note that the semantics for \code{which_multiple} and \code{count_multiple} is
#' different. \code{which_multiple} gives \code{TRUE} for all occurences of a
#' multiple edge except for one. Ie. if there are three \code{i-j} edges in the
#' graph then \code{which_multiple} returns \code{TRUE} for only two of them while
#' \code{count_multiple} returns \sQuote{3} for all three.
#' 
#' See the examples for getting rid of multiple edges while keeping their
#' original multiplicity as an edge attribute.
#' 
#' @aliases has.multiple is.loop is.multiple count.multiple count_multiple
#'   any_multiple which_loop
#' @param graph The input graph.
#' @param eids The edges to which the query is restricted. By default this is
#' all edges in the graph.
#' @return \code{any_multiple} returns a logical scalar.  \code{which_loop} and
#' \code{which_multiple} return a logical vector. \code{count_multiple} returns a
#' numeric vector.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{simplify}} to eliminate loop and multiple edges.
#' @export
#' @keywords graphs
#' @examples
#' 
#' # Loops
#' g <- graph( c(1,1,2,2,3,3,4,5) )
#' which_loop(g)
#' 
#' # Multiple edges
#' g <- barabasi.game(10, m=3, algorithm="bag")
#' any_multiple(g)
#' which_multiple(g)
#' count_multiple(g)
#' which_multiple(simplify(g))
#' all(count_multiple(simplify(g)) == 1)
#' 
#' # Direction of the edge is important
#' which_multiple(graph( c(1,2, 2,1) ))
#' which_multiple(graph( c(1,2, 2,1), dir=FALSE ))
#' 
#' # Remove multiple edges but keep multiplicity
#' g <- barabasi.game(10, m=3, algorithm="bag")
#' E(g)$weight <- count_multiple(g)
#' g <- simplify(g)
#' any(which_multiple(g))
#' E(g)$weight
#' 
which_multiple <- function(graph, eids=E(graph)) {

  if (!is_igraph(graph)) {
    stop("Not a graph object");
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_is_multiple", graph, as.igraph.es(graph, eids)-1,
        PACKAGE="igraph")
}

#' @export

count_multiple <- function(graph, eids=E(graph)) {

  if (!is_igraph(graph)) {
    stop("Not a graph object");
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_count_multiple", graph, as.igraph.es(graph, eids)-1,
        PACKAGE="igraph")
}



#' Breadth-first search
#' 
#' Breadth-first search is an algorithm to traverse a graph. We start from a
#' root vertex and spread along every edge \dQuote{simultaneously}.
#' 
#' 
#' The callback function must have the following arguments: \describe{
#' \item{graph}{The input graph is passed to the callback function here.}
#' \item{data}{A named numeric vector, with the following entries:
#' \sQuote{vid}, the vertex that was just visited, \sQuote{pred}, its
#' predecessor, \sQuote{succ}, its successor, \sQuote{rank}, the rank of the
#' current vertex, \sQuote{dist}, its distance from the root of the search
#' tree.} \item{extra}{The extra argument.} } See examples below on how to use
#' the callback function.
#'
#' @aliases graph.bfs
#' @param graph The input graph.
#' @param root Numeric vector, usually of length one. The root vertex, or root
#' vertices to start the search from.
#' @param neimode For directed graphs specifies the type of edges to follow.
#' \sQuote{out} follows outgoing, \sQuote{in} incoming edges. \sQuote{all}
#' ignores edge directions completely. \sQuote{total} is a synonym for
#' \sQuote{all}. This argument is ignored for undirected graphs.
#' @param unreachable Logical scalar, whether the search should visit the
#' vertices that are unreachable from the given root vertex (or vertices). If
#' \code{TRUE}, then additional searches are performed until all vertices are
#' visited.
#' @param restricted \code{NULL} (=no restriction), or a vector of vertices
#' (ids or symbolic names). In the latter case, the search is restricted to the
#' given vertices.
#' @param order Logical scalar, whether to return the ordering of the vertices.
#' @param rank Logical scalar, whether to return the rank of the vertices.
#' @param father Logical scalar, whether to return the father of the vertices.
#' @param pred Logical scalar, whether to return the predecessors of the
#' vertices.
#' @param succ Logical scalar, whether to return the successors of the
#' vertices.
#' @param dist Logical scalar, whether to return the distance from the root of
#' the search tree.
#' @param callback If not \code{NULL}, then it must be callback function. This
#' is called whenever a vertex is visited. See details below.
#' @param extra Additional argument to supply to the callback function.
#' @param rho The environment in which the callback function is evaluated.
#' @return A named list with the following entries: \item{root}{Numeric scalar.
#' The root vertex that was used as the starting point of the search.}
#' \item{neimode}{Character scalar. The \code{neimode} argument of the function
#' call. Note that for undirected graphs this is always \sQuote{all},
#' irrespectively of the supplied value.} \item{order}{Numeric vector. The
#' vertex ids, in the order in which they were visited by the search.}
#' \item{rank}{Numeric vector. The rank for each vertex.} \item{father}{Numeric
#' vector. The father of each vertex, i.e. the vertex it was discovered from.}
#' \item{pred}{Numeric vector. The previously visited vertex for each vertex,
#' or 0 if there was no such vertex.} \item{succ}{Numeric vector. The next
#' vertex that was visited after the current one, or 0 if there was no such
#' vertex.} \item{dist}{Numeric vector, for each vertex its distance from the
#' root of the search tree.}
#' 
#' Note that \code{order}, \code{rank}, \code{father}, \code{pred}, \code{succ}
#' and \code{dist} might be \code{NULL} if their corresponding argument is
#' \code{FALSE}, i.e. if their calculation is not requested.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{dfs}} for depth-first search.
#' @export
#' @keywords graphs
#' @examples
#' 
#' ## Two rings
#' bfs(make_ring(10) %du% make_ring(10), root=1, "out",
#'           order=TRUE, rank=TRUE, father=TRUE, pred=TRUE,
#'           succ=TRUE, dist=TRUE)
#' 
#' ## How to use a callback
#' f <- function(graph, data, extra) {
#'   print(data)
#'   FALSE
#' }
#' tmp <- bfs(make_ring(10) %du% make_ring(10), root=1, "out",
#'                  callback=f)
#' 
#' ## How to use a callback to stop the search
#' ## We stop after visiting all vertices in the initial component
#' f <- function(graph, data, extra) {
#'  data['succ'] == -1
#' }
#' bfs(make_ring(10) %du% make_ring(10), root=1, callback=f)
#' 
#' 
bfs <- function(graph, root, neimode=c("out", "in", "all", "total"),
                      unreachable=TRUE, restricted=NULL,
                      order=TRUE, rank=FALSE, father=FALSE,
                      pred=FALSE, succ=FALSE, dist=FALSE,
                      callback=NULL, extra=NULL, rho=parent.frame()) {

  if (!is_igraph(graph)) {
    stop("Not a graph object");
  }

  if (length(root)==1) {
    root <- as.igraph.vs(graph, root)-1
    roots <- NULL    
  } else {
    roots <- as.igraph.vs(graph, root)-1
    root <- 0      # ignored anyway
  }
  neimode <- switch(igraph.match.arg(neimode),
                    "out"=1, "in"=2, "all"=3, "total"=3)
  unreachable <- as.logical(unreachable)
  if (!is.null(restricted)) { restricted <- as.igraph.vs(graph, restricted) }
  if (!is.null(callback)) { callback <- as.function(callback) }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_bfs", graph, root, roots, neimode, unreachable,
               restricted,
               as.logical(order), as.logical(rank), as.logical(father),
               as.logical(pred), as.logical(succ), as.logical(dist),
               callback, extra, rho,
               PACKAGE="igraph")
  
  if (order)  res$order  <- res$order+1
  if (rank)   res$rank   <- res$rank+1
  if (father) res$father <- res$father+1
  if (pred)   res$pred   <- res$pred+1
  if (succ)   res$succ   <- res$succ+1
  res
}



#' Depth-first search
#' 
#' Depth-first search is an algorithm to traverse a graph. It starts from a
#' root vertex and tries to go quickly as far from as possible.
#' 
#' The callback functions must have the following arguments: \describe{
#' \item{graph}{The input graph is passed to the callback function here.}
#' \item{data}{A named numeric vector, with the following entries:
#' \sQuote{vid}, the vertex that was just visited and \sQuote{dist}, its
#' distance from the root of the search tree.} \item{extra}{The extra
#' argument.} } See examples below on how to use the callback functions.
#'
#' @aliases graph.dfs
#' @param graph The input graph.
#' @param root The single root vertex to start the search from.
#' @param neimode For directed graphs specifies the type of edges to follow.
#' \sQuote{out} follows outgoing, \sQuote{in} incoming edges. \sQuote{all}
#' ignores edge directions completely. \sQuote{total} is a synonym for
#' \sQuote{all}. This argument is ignored for undirected graphs.
#' @param unreachable Logical scalar, whether the search should visit the
#' vertices that are unreachable from the given root vertex (or vertices). If
#' \code{TRUE}, then additional searches are performed until all vertices are
#' visited.
#' @param order Logical scalar, whether to return the DFS ordering of the
#' vertices.
#' @param order.out Logical scalar, whether to return the ordering based on
#' leaving the subtree of the vertex.
#' @param father Logical scalar, whether to return the father of the vertices.
#' @param dist Logical scalar, whether to return the distance from the root of
#' the search tree.
#' @param in.callback If not \code{NULL}, then it must be callback function.
#' This is called whenever a vertex is visited. See details below.
#' @param out.callback If not \code{NULL}, then it must be callback function.
#' This is called whenever the subtree of a vertex is completed by the
#' algorithm. See details below.
#' @param extra Additional argument to supply to the callback function.
#' @param rho The environment in which the callback function is evaluated.
#' @return A named list with the following entries: \item{root}{Numeric scalar.
#' The root vertex that was used as the starting point of the search.}
#' \item{neimode}{Character scalar. The \code{neimode} argument of the function
#' call. Note that for undirected graphs this is always \sQuote{all},
#' irrespectively of the supplied value.} \item{order}{Numeric vector. The
#' vertex ids, in the order in which they were visited by the search.}
#' \item{order.out}{Numeric vector, the vertex ids, in the order of the
#' completion of their subtree.} \item{father}{Numeric vector. The father of
#' each vertex, i.e. the vertex it was discovered from.} \item{dist}{Numeric
#' vector, for each vertex its distance from the root of the search tree.}
#' 
#' Note that \code{order}, \code{order.out}, \code{father}, and \code{dist}
#' might be \code{NULL} if their corresponding argument is \code{FALSE}, i.e.
#' if their calculation is not requested.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{bfs}} for breadth-first search.
#' @export
#' @keywords graphs
#' @examples
#' 
#' ## A graph with two separate trees
#' dfs(make_tree(10) %du% make_tree(10), root=1, "out",
#'           TRUE, TRUE, TRUE, TRUE)
#' 
#' ## How to use a callback
#' f.in <- function(graph, data, extra) {
#'   cat("in:", paste(collapse=", ", data), "\n")
#'   FALSE
#' }
#' f.out <- function(graph, data, extra) {
#'   cat("out:", paste(collapse=", ", data), "\n")
#'   FALSE
#' }
#' tmp <- dfs(make_tree(10), root=1, "out",
#'                  in.callback=f.in, out.callback=f.out)
#' 
#' ## Terminate after the first component, using a callback
#' f.out <- function(graph, data, extra) {
#'  data['vid'] == 1
#' }
#' tmp <- dfs(make_tree(10) %du% make_tree(10), root=1,
#'                  out.callback=f.out)
#' 
#' 
dfs <- function(graph, root, neimode=c("out", "in", "all", "total"),
                      unreachable=TRUE,
                      order=TRUE, order.out=FALSE, father=FALSE, dist=FALSE,
                      in.callback=NULL, out.callback=NULL, extra=NULL,
                      rho=parent.frame()) {

  if (!is_igraph(graph)) {
    stop("Not a graph object");
  }

  root <- as.igraph.vs(graph, root)-1
  neimode <- switch(igraph.match.arg(neimode),
                    "out"=1, "in"=2, "all"=3, "total"=3)
  unreachable <- as.logical(unreachable)
  if (!is.null(in.callback)) { in.callback <- as.function(in.callback) }
  if (!is.null(out.callback)) { out.callback <- as.function(out.callback) }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_dfs", graph, root, neimode, unreachable,
               as.logical(order), as.logical(order.out), as.logical(father),
               as.logical(dist), in.callback, out.callback, extra, rho,
               PACKAGE="igraph")
  
  if (order)     res$order     <- res$order+1
  if (order.out) res$order.out <- res$order.out+1
  if (father)    res$father    <- res$father+1
  res
}

#' @rdname betweenness
#' @param e The edges for which the edge betweenness will be calculated.
#' @export

edge_betweenness <- function(graph, e=E(graph),
                             directed=TRUE, weights=NULL) {
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  e <- as.igraph.es(graph, e)
  directed <- as.logical(directed)
  if (is.null(weights) && "weight" %in% edge_attr_names(graph)) { 
  weights <- E(graph)$weight 
  } 
  if (!is.null(weights) && any(!is.na(weights))) { 
  weights <- as.numeric(weights) 
  } else { 
  weights <- NULL 
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_edge_betweenness", graph, directed, weights,
        PACKAGE="igraph")
  res[as.numeric(e)]
}

#' @export

estimate_edge_betweenness <- function(graph, e=E(graph),
                                      directed=TRUE, cutoff, weights=NULL) {
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  e <- as.igraph.es(graph, e)
  directed <- as.logical(directed)
  cutoff <- as.numeric(cutoff)
  if (is.null(weights) && "weight" %in% edge_attr_names(graph)) { 
  weights <- E(graph)$weight 
  } 
  if (!is.null(weights) && any(!is.na(weights))) { 
  weights <- as.numeric(weights) 
  } else { 
  weights <- NULL 
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_edge_betweenness_estimate", graph, directed, cutoff, weights,
        PACKAGE="igraph")
  res[as.numeric(e)]
}



#' Connected components of a graph
#' 
#' Calculate the maximal (weakly or strongly) connected components of a graph
#' 
#' \code{is_connected} decides whether the graph is weakly or strongly
#' connected.
#' 
#' \code{components} finds the maximal (weakly or strongly) connected components
#' of a graph.
#' 
#' \code{count_components} does almost the same as \code{components} but returns only
#' the number of clusters found instead of returning the actual clusters.
#' 
#' \code{component_distribution} creates a histogram for the maximal connected
#' component sizes.
#' 
#' The weakly connected components are found by a simple breadth-first search.
#' The strongly connected components are implemented by two consecutive
#' depth-first searches.
#' 
#' @aliases no.clusters clusters is.connected cluster.distribution components
#'   count_components is_connected
#' @param graph The graph to analyze.
#' @param mode Character string, either \dQuote{weak} or \dQuote{strong}.  For
#' directed graphs \dQuote{weak} implies weakly, \dQuote{strong} strongly
#' connected components to search. It is ignored for undirected graphs.
#' @param \dots Additional attributes to pass to \code{cluster}, right now only
#' \code{mode} makes sense.
#' @return For \code{is_connected} a logical constant.
#' 
#' For \code{components} a named list with three components:
#' \item{membership}{numeric vector giving the cluster id to which each vertex
#' belongs.} \item{csize}{numeric vector giving the sizes of the clusters.}
#' \item{no}{numeric constant, the number of clusters.}
#' 
#' For \code{count_components} an integer constant is returned.
#' 
#' For \code{component_distribution} a numeric vector with the relative
#' frequencies. The length of the vector is the size of the largest component
#' plus one. Note that (for currently unknown reasons) the first element of the
#' vector is the number of clusters of size zero, so this is always zero.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{subcomponent}}, \code{\link{groups}}
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- sample_gnp(20, 1/20)
#' clu <- components(g)
#' groups(clu)
#' 
components <- function(graph, mode=c("weak", "strong")) {
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  mode <- switch(igraph.match.arg(mode), "weak"=1, "strong"=2)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_clusters", graph, mode,
        PACKAGE="igraph")
  res$membership <- res$membership + 1
  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    names(res$membership) <- V(graph)$name
  }
  res
}



#' Convert a general graph into a forest
#' 
#' Perform a breadth-first search on a graph and convert it into a tree or
#' forest by replicating vertices that were found more than once.
#' 
#' A forest is a graph, whose components are trees.
#' 
#' The \code{roots} vector can be calculated by simply doing a topological sort
#' in all components of the graph, see the examples below.
#' 
#' @aliases unfold.tree
#' @param graph The input graph, it can be either directed or undirected.
#' @param mode Character string, defined the types of the paths used for the
#' breadth-first search. \dQuote{out} follows the outgoing, \dQuote{in} the
#' incoming edges, \dQuote{all} and \dQuote{total} both of them. This argument
#' is ignored for undirected graphs.
#' @param roots A vector giving the vertices from which the breadth-first
#' search is performed. Typically it contains one vertex per component.
#' @return A list with two components: \item{tree}{The result, an \code{igraph}
#' object, a tree or a forest.} \item{vertex_index}{A numeric vector, it gives
#' a mapping from the vertices of the new graph to the vertices of the old
#' graph.}
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- make_tree(10) %du% make_tree(10)
#' V(g)$id <- seq_len(vcount(g))-1
#' roots <- sapply(decompose(g), function(x) {
#'             V(x)$id[ topo_sort(x)[1]+1 ] })
#' tree <- unfold_tree(g, roots=roots)
#' 
unfold_tree <- function(graph, mode=c("all", "out", "in", "total"), roots) {
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  mode <- switch(igraph.match.arg(mode), "out"=1, "in"=2, "all"=3, "total"=3)
  roots <- as.igraph.vs(graph, roots)-1

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_unfold_tree", graph, mode, roots,
        PACKAGE="igraph")
  res
}



#' Closeness centrality of vertices
#' 
#' Cloness centrality measures how many steps is required to access every other
#' vertex from a given vertex.
#' 
#' The closeness centrality of a vertex is defined by the inverse of the
#' average length of the shortest paths to/from all the other vertices in the
#' graph:
#' 
#' \deqn{\frac{1}{\sum_{i\ne v} d_vi}}{1/sum( d(v,i), i != v)}
#' 
#' If there is no (directed) path between vertex \eqn{v}{\code{v}} and
#' \eqn{i}{\code{i}} then the total number of vertices is used in the formula
#' instead of the path length.
#' 
#' \code{estimate_closeness} only considers paths of length \code{cutoff} or
#' smaller, this can be run for larger graphs, as the running time is not
#' quadratic (if \code{cutoff} is small). If \code{cutoff} is zero or negative
#' then the function calculates the exact closeness scores.
#' 
#' @aliases closeness closeness.estimate estimate_closeness
#' @param graph The graph to analyze.
#' @param vids The vertices for which closeness will be calculated.
#' @param mode Character string, defined the types of the paths used for
#' measuring the distance in directed graphs. \dQuote{in} measures the paths
#' \emph{to} a vertex, \dQuote{out} measures paths \emph{from} a vertex,
#' \emph{all} uses undirected paths. This argument is ignored for undirected
#' graphs.
#' @param normalized Logical scalar, whether to calculate the normalized
#' closeness. Normalization is performed by multiplying the raw closeness by
#' \eqn{n-1}, where \eqn{n} is the number of vertices in the graph.
#' @param weights Optional positive weight vector for calculating weighted
#' closeness. If the graph has a \code{weight} edge attribute, then this is
#' used by default.
#' @return Numeric vector with the closeness values of all the vertices in
#' \code{v}.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{betweenness}}, \code{\link{degree}}
#' @references Freeman, L.C. (1979). Centrality in Social Networks I:
#' Conceptual Clarification. \emph{Social Networks}, 1, 215-239.
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- make_ring(10)
#' g2 <- make_star(10)
#' closeness(g)
#' closeness(g2, mode="in")
#' closeness(g2, mode="out")
#' closeness(g2, mode="all")
#' 
closeness <- function(graph, vids=V(graph),
                      mode=c("out", "in", "all", "total"), weights=NULL,
                      normalized=FALSE) {
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  vids <- as.igraph.vs(graph, vids)
  mode <- switch(igraph.match.arg(mode), "out"=1, "in"=2, "all"=3, "total"=3)
  if (is.null(weights) && "weight" %in% edge_attr_names(graph)) { 
  weights <- E(graph)$weight 
  } 
  if (!is.null(weights) && any(!is.na(weights))) { 
  weights <- as.numeric(weights) 
  } else { 
  weights <- NULL 
  }
  normalized <- as.logical(normalized)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_closeness", graph, vids-1, mode, weights,
               normalized, PACKAGE="igraph")
  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    names(res) <- V(graph)$name[vids]
  }
  res
}

#' @rdname closeness
#' @param cutoff The maximum path length to consider when calculating the
#' betweenness. If zero or negative then there is no such limit.
#' @export

estimate_closeness <- function(graph, vids=V(graph), mode=c("out", "in", "all", "total"), cutoff, weights=NULL, normalized=FALSE) {
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  vids <- as.igraph.vs(graph, vids)
  mode <- switch(igraph.match.arg(mode), "out"=1, "in"=2, "all"=3, "total"=3)
  cutoff <- as.numeric(cutoff)
  if (is.null(weights) && "weight" %in% edge_attr_names(graph)) { 
  weights <- E(graph)$weight 
  } 
  if (!is.null(weights) && any(!is.na(weights))) { 
  weights <- as.numeric(weights) 
  } else { 
  weights <- NULL 
  }
  normalized <- as.logical(normalized)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_closeness_estimate", graph, vids-1, mode, cutoff, weights, normalized,
        PACKAGE="igraph")
  if (igraph_opt("add.vertex.names") && is_named(graph)) { 
  names(res) <- vertex_attr(graph, "name", vids) 
  }
  res
}


#' Graph Laplacian
#' 
#' The Laplacian of a graph.
#' 
#' The Laplacian Matrix of a graph is a symmetric matrix having the same number
#' of rows and columns as the number of vertices in the graph and element (i,j)
#' is d[i], the degree of vertex i if if i==j, -1 if i!=j and there is an edge
#' between vertices i and j and 0 otherwise.
#' 
#' A normalized version of the Laplacian Matrix is similar: element (i,j) is 1
#' if i==j, -1/sqrt(d[i] d[j]) if i!=j and there is an edge between vertices i
#' and j and 0 otherwise.
#' 
#' The weighted version of the Laplacian simply works with the weighted degree
#' instead of the plain degree. I.e. (i,j) is d[i], the weighted degree of
#' vertex i if if i==j, -w if i!=j and there is an edge between vertices i and
#' j with weight w, and 0 otherwise. The weighted degree of a vertex is the sum
#' of the weights of its adjacent edges.
#' 
#' @param graph The input graph.
#' @param normalized Whether to calculate the normalized Laplacian. See
#' definitions below.
#' @param weights An optional vector giving edge weights for weighted Laplacian
#' matrix. If this is \code{NULL} and the graph has an edge attribute called
#' \code{weight}, then it will be used automatically. Set this to \code{NA} if
#' you want the unweighted Laplacian on a graph that has a \code{weight} edge
#' attribute.
#' @param sparse Logical scalar, whether to return the result as a sparse
#' matrix. The \code{Matrix} package is required for sparse matrices.
#' @return A numeric matrix.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- make_ring(10)
#' laplacian_matrix(g)
#' laplacian_matrix(g, norm=TRUE)
#' laplacian_matrix(g, norm=TRUE, sparse=FALSE)
#' 
laplacian_matrix <- function(graph, normalized=FALSE, weights=NULL,
                            sparse=igraph_opt("sparsematrices")) {
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  normalized <- as.logical(normalized)
  if (is.null(weights) && "weight" %in% edge_attr_names(graph)) { 
    weights <- E(graph)$weight 
  } 
  if (!is.null(weights) && any(!is.na(weights))) { 
    weights <- as.numeric(weights) 
  } else { 
    weights <- NULL 
  }
  sparse <- as.logical(sparse)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_laplacian", graph, normalized, weights, sparse,
               PACKAGE="igraph")
  if (sparse) {
    res <- igraph.i.spMatrix(res)
  }
  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    rownames(res) <- colnames(res) <- V(graph)$name
  }
  res
}

#' @keywords graphs
 
is_matching <- function(graph, matching, types=NULL) {
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  if (is.null(types) && "type" %in% vertex_attr_names(graph)) { 
    types <- V(graph)$type 
  } 
  if (!is.null(types)) { 
    types <- as.logical(types) 
  }
  matching <- as.igraph.vs(graph, matching, na.ok=TRUE)-1
  matching[ is.na(matching) ] <- -1

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_is_matching", graph, types, matching,
        PACKAGE="igraph")

  res
}

#' @keywords graphs
 
is_max_matching <- function(graph, matching, types=NULL) {
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  if (is.null(types) && "type" %in% vertex_attr_names(graph)) { 
    types <- V(graph)$type 
  } 
  if (!is.null(types)) { 
    types <- as.logical(types) 
  }
  matching <- as.igraph.vs(graph, matching, na.ok=TRUE)-1
  matching[ is.na(matching) ] <- -1

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_is_maximal_matching", graph, types, matching,
        PACKAGE="igraph")

  res
}

#' @keywords graphs
 
max_bipartite_match <- function(graph, types=NULL, weights=NULL,
                                       eps=.Machine$double.eps) {
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  if (is.null(types) && "type" %in% vertex_attr_names(graph)) { 
    types <- V(graph)$type 
  } 
  if (!is.null(types)) { 
    types <- as.logical(types) 
  }
  if (is.null(weights) && "weight" %in% edge_attr_names(graph)) { 
    weights <- E(graph)$weight 
  } 
  if (!is.null(weights) && any(!is.na(weights))) { 
    weights <- as.numeric(weights) 
  } else { 
    weights <- NULL 
  }
  eps <- as.numeric(eps)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_maximum_bipartite_matching", graph, types, weights,
               eps,
               PACKAGE="igraph")

  res$matching[ res$matching==0 ] <- NA
  if (igraph_opt("add.vertex.names") && is_named(graph)) {
    res$matching <- V(graph)$name[res$matching]
    names(res$matching) <- V(graph)$name
  }
  res
}
