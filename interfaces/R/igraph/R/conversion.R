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

get.adjacency.dense <- function(graph, type=c("both", "upper", "lower"),
                                attr=NULL, edges=FALSE, names=TRUE) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  
  type <- igraph.match.arg(type)
  type <- switch(type, "upper"=0, "lower"=1, "both"=2)
  
  if (edges || is.null(attr)) {    
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    res <- .Call("R_igraph_get_adjacency", graph, as.numeric(type),
                 as.logical(edges), PACKAGE="igraph")
  } else {
    attr <- as.character(attr)
    if (! attr %in% edge_attr_names(graph)) {
      stop("no such edge attribute")
    }
    exattr <- edge_attr(graph, attr)
    if (is.logical(exattr)) {
      res <- matrix(FALSE, nrow=vcount(graph), ncol=vcount(graph))
    } else if (is.character(exattr)) {
      res <- matrix("", nrow=vcount(graph), ncol=vcount(graph))
    } else if (is.numeric(exattr)) {
      res <- matrix(0, nrow=vcount(graph), ncol=vcount(graph))
    } else {
      stop("Sparse matrices must be either numeric or logical,",
           "and the edge attribute is not")
    }
    if (is_directed(graph)) {
      for (i in seq(length=ecount(graph))) {
        e <- get.edge(graph, i)
        res[ e[1], e[2] ] <- edge_attr(graph, attr, i)
      }
    } else {
      if (type==0) {
        ## upper
        for (i in seq(length=ecount(graph))) {
          e <- get.edge(graph, i)
          res[ min(e), max(e) ] <- edge_attr(graph, attr, i)
        }        
      } else if (type==1) {
        ## lower
        for (i in seq(length=ecount(graph))) {
          e <- get.edge(graph, i)
          res[ max(e), min(e) ] <- edge_attr(graph, attr, i)
        }        
      } else if (type==2) {
        ## both
        for (i in seq(length=ecount(graph))) {
          e <- get.edge(graph, i)
          res[ e[1], e[2] ] <- edge_attr(graph, attr, i)
          if (e[1] != e[2]) {
            res[ e[2], e[1] ] <- edge_attr(graph, attr, i)
          }
        }
      }
    }
  }

  if (names && "name" %in% vertex_attr_names(graph)) {
    colnames(res) <- rownames(res) <- V(graph)$name
  }
  
  res  
}

get.adjacency.sparse <- function(graph, type=c("both", "upper", "lower"),
                                 attr=NULL, edges=FALSE, names=TRUE) {

  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  type <- igraph.match.arg(type)

  vc <- vcount(graph)
  
  el <- as_edgelist(graph, names=FALSE)
  if (edges) {
    value <- seq_len(nrow(el))
  } else if (!is.null(attr)) {
    attr <- as.character(attr)
    if (!attr %in% edge_attr_names(graph)) {
      stop("no such edge attribute")
    }
    value <- edge_attr(graph, name=attr)
    if (!is.numeric(value) && !is.logical(value)) {
      stop("Sparse matrices must be either numeric or logical,",
           "and the edge attribute is not")
    }
  } else {
    value <- rep(1, nrow(el))
  }

  if (is_directed(graph)) {
    res <- Matrix::sparseMatrix(dims=c(vc, vc), i=el[,1], j=el[,2], x=value)
  } else {
    if (type=="upper") {
      ## upper
      res <- Matrix::sparseMatrix(dims=c(vc, vc), i=pmin(el[,1],el[,2]),
                          j=pmax(el[,1],el[,2]), x=value)
    } else if (type=="lower") {
      ## lower
      res <- Matrix::sparseMatrix(dims=c(vc, vc), i=pmax(el[,1],el[,2]),
                          j=pmin(el[,1],el[,2]), x=value)
    } else if (type=="both") {
      ## both
      res <- Matrix::sparseMatrix(dims=c(vc, vc), i=pmin(el[,1],el[,2]),
                          j=pmax(el[,1],el[,2]), x=value, symmetric=TRUE)
      res <- as(res, "dgCMatrix")
    }
  }

  if (names && "name" %in% vertex_attr_names(graph)) {
    colnames(res) <- rownames(res) <- V(graph)$name
  }

  res
}

#' @export

as_adj <- function(graph, type=c("both", "upper", "lower"),
                          attr=NULL, edges=FALSE, names=TRUE, 
                          sparse=igraph_opt("sparsematrices")) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  if (!sparse) {
    get.adjacency.dense(graph, type=type, attr=attr, edges=edges, names=names)
  } else {
    get.adjacency.sparse(graph, type=type, attr=attr, edges=edges, names=names)
  }  
}

#' @export

as_adjacency_matrix <- as_adj

#' @export

as_edgelist <- function(graph, names=TRUE) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- matrix(.Call("R_igraph_get_edgelist", graph, TRUE,
                      PACKAGE="igraph"), ncol=2)
  res <- res+1
  if (names && "name" %in% vertex_attr_names(graph)) {
    res <- matrix(V(graph)$name[ res ], ncol=2)
  }

  res
}



#' Convert between directed and undirected graphs
#' 
#' \code{as.directed} converts an undirected graph to directed,
#' \code{as.undirected} does the opposite, it converts a directed graph to
#' undirected.
#' 
#' Conversion algorithms for \code{as.directed}: \describe{
#' \item{list("arbitrary")}{The number of edges in the graph stays the same, an
#' arbitrarily directed edge is created for each undirected edge.}
#' \item{list("mutual")}{Two directed edges are created for each undirected
#' edge, one in each direction.} }
#' 
#' Conversion algorithms for \code{as.undirected}: \describe{
#' \item{list("each")}{The number of edges remains constant, an undirected edge
#' is created for each directed one, this version might create graphs with
#' multiple edges.} \item{list("collapse")}{One undirected edge will be created
#' for each pair of vertices which are connected with at least one directed
#' edge, no multiple edges will be created.} \item{list("mutual")}{One
#' undirected edge will be created for each pair of mutual edges. Non-mutual
#' edges are ignored. This mode might create multiple edges if there are more
#' than one mutual edge pairs between the same pair of vertices.  } }
#' 
#' @aliases as.directed as.undirected
#' @param graph The graph to convert.
#' @param mode Character constant, defines the conversion algorithm. For
#' \code{as.directed} it can be \code{mutual} or \code{arbitrary}. For
#' \code{as.undirected} it can be \code{each}, \code{collapse} or
#' \code{mutual}. See details below.
#' @return A new graph object.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{simplify}} for removing multiple and/or loop edges from
#' a graph.
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- ring(10)
#' as.directed(g, "mutual")
#' g2 <- star(10)
#' as.undirected(g)
#' 
#' # Combining edge attributes
#' g3 <- ring(10, directed=TRUE, mutual=TRUE)
#' E(g3)$weight <- seq_len(ecount(g3))
#' ug3 <- as.undirected(g3)
#' print(ug3, e=TRUE)
#' \dontrun{
#'   x11(width=10, height=5)
#'   layout(rbind(1:2))
#'   plot( g3, layout=layout_in_circle, edge.label=E(g3)$weight)
#'   plot(ug3, layout=layout_in_circle, edge.label=E(ug3)$weight)
#' }
#' 
#' g4 <- graph(c(1,2, 3,2,3,4,3,4, 5,4,5,4,
#'               6,7, 7,6,7,8,7,8, 8,7,8,9,8,9,
#'               9,8,9,8,9,9, 10,10,10,10))
#' E(g4)$weight <- seq_len(ecount(g4))
#' ug4 <- as.undirected(g4, mode="mutual",
#'               edge.attr.comb=list(weight=length))
#' print(ug4, e=TRUE)
#' 
as.directed <- function(graph, mode=c("mutual", "arbitrary")) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "arbitrary"=0, "mutual"=1)
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_to_directed", graph, as.numeric(mode),
        PACKAGE="igraph")
}

#' @rdname as.directed
#' @param edge.attr.comb Specifies what to do with edge attributes, if
#' \code{mode="collapse"} or \code{mode="mutual"}.  In these cases many edges
#' might be mapped to a single one in the new graph, and their attributes are
#' combined. Please see \code{\link{attribute.combination}} for details on
#' this.
#' @export

as.undirected <- function(graph, mode=c("collapse", "each", "mutual"), edge.attr.comb=igraph_opt("edge.attr.comb")) {
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  mode <- switch(igraph.match.arg(mode), "collapse"=1, "each"=0, "mutual"=2)
  edge.attr.comb <- igraph.i.attribute.combination(edge.attr.comb)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_to_undirected", graph, mode, edge.attr.comb,
        PACKAGE="igraph")

  res
}


#' Adjacency lists
#' 
#' Create adjacency lists from a graph, either for adjacent edges or for
#' neighboring vertices
#' 
#' \code{as_adj_list} returns a list of numeric vectors, which include the ids
#' of neighbor vertices (according to the \code{mode} argument) of all
#' vertices.
#' 
#' \code{as_adj_edge_list} returns a list of numeric vectors, which include the
#' ids of adjacent edgs (according to the \code{mode} argument) of all
#' vertices.
#' 
#' @aliases as_adj_list get.adjedgelist
#' @param graph The input graph.
#' @param mode Character scalar, it gives what kind of adjacent edges/vertices
#' to include in the lists. \sQuote{\code{out}} is for outgoing edges/vertices,
#' \sQuote{\code{in}} is for incoming edges/vertices, \sQuote{\code{all}} is
#' for both. This argument is ignored for undirected graphs.
#' @return A list of numeric vectors.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{as_edgelist}}, \code{\link{as_adj}}
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- ring(10)
#' as_adj_list(g)
#' as_adj_edge_list(g)
#' 
as_adj_list <- function(graph, mode=c("all", "out", "in", "total")) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  mode <- igraph.match.arg(mode)
  mode <- as.numeric(switch(mode, "out"=1, "in"=2, "all"=3, "total"=3))
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_get_adjlist", graph, mode,
               PACKAGE="igraph")
  res <- lapply(res, function(x) x+1)
  if (is_named(graph)) names(res) <- V(graph)$name
  res
}

#' @rdname as_adj_list
#' @aliases get.adjlist
#' @export

as_adj_edge_list <- function(graph, mode=c("all", "out", "in", "total")) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  mode <- igraph.match.arg(mode)
  mode <- as.numeric(switch(mode, "out"=1, "in"=2, "all"=3, "total"=3))
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_get_adjedgelist", graph, mode,
               PACKAGE="igraph")
  res <- lapply(res, function(x) x+1)
  if (is_named(graph)) names(res) <- V(graph)$name
  res
}

#' @export

graph_from_graphnel <- function(graphNEL, name=TRUE, weight=TRUE,
                                 unlist.attrs=TRUE) {

  if (! "graph" %in% .packages()) {
    library(graph, pos="package:base")
  }

  if (!inherits(graphNEL, "graphNEL")) {
    stop("Not a graphNEL graph")
  }
  
  al <- lapply(edgeL(graphNEL), "[[", "edges")
  if (edgemode(graphNEL)=="undirected") {
    al <- mapply(SIMPLIFY=FALSE, seq_along(al), al, FUN=function(n, l) {
      c(l, rep(n, sum(l==n)))
    })
  }
  mode <- if (edgemode(graphNEL)=="directed") "out" else "all"
  g <- graph_from_adj_list(al, mode=mode, duplicate=TRUE)
  if (name) {
    V(g)$name <- nodes(graphNEL)
  }

  ## Graph attributes
  g.n <- names(graphNEL@graphData)
  g.n <- g.n [ g.n != "edgemode" ]
  for (n in g.n) {
    g <- set_graph_attr(g, n, graphNEL@graphData[[n]])
  }
  
  ## Vertex attributes
  v.n <- names(nodeDataDefaults(graphNEL))
  for (n in v.n) {
    val <- unname(nodeData(graphNEL, attr=n))
    if (unlist.attrs && all(sapply(val, length)==1)) { val <- unlist(val) }
    g <- set_vertex_attr(g, n, value=val)
  }

  ## Edge attributes
  e.n <- names(edgeDataDefaults(graphNEL))
  if (!weight) { e.n <- e.n [ e.n != "weight" ] }
  if (length(e.n) > 0) {
    el <- as_edgelist(g)
    el <- paste(sep="|", el[,1], el[,2])
    for (n in e.n) {
      val <- unname(edgeData(graphNEL, attr=n)[el])
      if (unlist.attrs && all(sapply(val, length)==1)) { val <- unlist(val) }
      g <- set_edge_attr(g, n, value=val)
    }
  }
  
  g 
}

#' @export

as_graphnel <- function(graph) {

  if (!is_igraph(graph)) {
    stop("Not an igraph graph")
  }
  
  if (! "graph" %in% .packages()) {
    library(graph, pos="package:base")
  }

  if ("name" %in% vertex_attr_names(graph) &&
      is.character(V(graph)$name)) {
    name <- V(graph)$name
  } else {
    name <- as.character(seq(vcount(graph)))    
  }

  edgemode <- if (is_directed(graph)) "directed" else "undirected"  

  if ("weight" %in% edge_attr_names(graph) &&
      is.numeric(E(graph)$weight)) {
    al <- as_adj_edge_list(graph, "out")
    for (i in seq(along=al)) {
      edges <- get.edges(graph, al[[i]])
      edges <- ifelse( edges[,2]==i, edges[,1], edges[,2])
      weights <- E(graph)$weight[al[[i]]]
      al[[i]] <- list(edges=edges, weights=weights)
    }
  } else {
    al <- as_adj_list(graph, "out")
    al <- lapply(al, function(x) list(edges=x))
  }  
  
  names(al) <- name
  res <- new("graphNEL", nodes=name, edgeL=al, edgemode=edgemode)

  ## Add graph attributes (other than 'directed')
  ## Are this "officially" supported at all?

  g.n <- graph_attr_names(graph)
  if ("directed" %in% g.n) {
    warning("Cannot add graph attribute `directed'")
    g.n <- g.n[ g.n != "directed" ]
  }
  for (n in g.n) {
    res@graphData[[n]] <- graph_attr(graph, n)
  }

  ## Add vertex attributes (other than 'name', that is already
  ## added as vertex names)
  
  v.n <- vertex_attr_names(graph)
  v.n <- v.n[ v.n != "name" ]
  for (n in v.n) {
    nodeDataDefaults(res, attr=n) <- NA
    nodeData(res, attr=n) <- vertex_attr(graph, n)
  }

  ## Add edge attributes (other than 'weight')
  
  e.n <- edge_attr_names(graph)
  e.n <- e.n[ e.n != "weight" ]
  if (length(e.n) > 0) {
    el <- as_edgelist(graph)
    el <- paste(sep="|", el[,1], el[,2])
    for (n in e.n) {
      edgeDataDefaults(res, attr=n) <- NA
      res@edgeData@data[el] <- mapply(function(x,y) {
        xx <- c(x,y); names(xx)[length(xx)] <- n; xx },
                                      res@edgeData@data[el],
                                      edge_attr(graph, n),
                                      SIMPLIFY=FALSE)
    }
  }
  
  res
}

get.incidence.dense <- function(graph, types, names, attr) {

  if (is.null(attr)) {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    ## Function call
    res <- .Call("R_igraph_get_incidence", graph, types,
                 PACKAGE="igraph")

    if (names && "name" %in% vertex_attr_names(graph)) {
      rownames(res$res) <- V(graph)$name[ res$row_ids+1 ]
      colnames(res$res) <- V(graph)$name[ res$col_ids+1 ]
    } else {
      rownames(res$res) <- res$row_ids+1
      colnames(res$res) <- res$col_ids+1
    }
    res$res
    
  } else {

    attr <- as.character(attr)
    if (!attr %in% edge_attr_names(graph)) {
      stop("no such edge attribute")
    }

    vc <- vcount(graph)
    n1 <- sum(!types)
    n2 <- vc-n1    
    res <- matrix(0, n1, n2)

    recode <- numeric(vc)
    recode[!types] <- seq_len(n1)
    recode[types]  <- seq_len(n2)
    
    for (i in seq(length=ecount(graph))) {
      eo <- get.edge(graph, i)
      e <- recode[eo]
      if (!types[eo[1]]) {
        res[ e[1], e[2] ] <- edge_attr(graph, attr, i)
      } else{
        res[ e[2], e[1] ] <- edge_attr(graph, attr, i)
      }
    }

    if (names && "name" %in% vertex_attr_names(graph)) {
      rownames(res) <- V(graph)$name[ which(!types) ]
      colnames(res) <- V(graph)$name[ which( types) ]
    } else {
      rownames(res) <- which(!types)
      colnames(res) <- which(types)
    }

    res
  }
}

get.incidence.sparse <- function(graph, types, names, attr) {

  vc <- vcount(graph)
  if (length(types) != vc) {
    stop("Invalid types vector")
  }

  el <- as_edgelist(graph, names=FALSE)
  if (any(types[el[,1]] == types[el[,2]])) {
    stop("Invalid types vector, not a bipartite graph")
  }

  n1 <- sum(!types)
  n2 <- vc-n1

  recode <- numeric(vc)
  recode[!types] <- seq_len(n1)
  recode[types]  <- seq_len(n2) + n1

  el[,1] <- recode[el[,1]]
  el[,2] <- recode[el[,2]]

  change <- el[,1] > n1
  el[change,] <- el[change,2:1]
  el[,2] <- el[,2]-n1

  if (!is.null(attr)) {
    attr <- as.character(attr)
    if (!attr %in% edge_attr_names(graph)) {
      stop("no such edge attribute")
    }
    value <- edge_attr(graph, name=attr)
  } else { 
    value <- rep(1, nrow(el))
  }

  res <- Matrix::spMatrix(n1, n2, i=el[,1], j=el[,2], x=value)

  if (names && "name" %in% vertex_attr_names(graph)) {
    rownames(res) <- V(graph)$name[which(!types)]
    colnames(res) <- V(graph)$name[which(types)]
  } else {
    rownames(res) <- which(!types)
    colnames(res) <- which(types)
  }
  res
}



#' Incidence matrix of a bipartite graph
#' 
#' This function can return a sparse or dense incidence matrix of a bipartite
#' network. The incidence matrix is an \eqn{n} times \eqn{m} matrix, \eqn{n}
#' and \eqn{m} are the number of vertices of the two kinds.
#' 
#' Bipartite graphs have a \code{type} vertex attribute in igraph, this is
#' boolean and \code{FALSE} for the vertices of the first kind and \code{TRUE}
#' for vertices of the second kind.
#' 
#' @param graph The input graph. The direction of the edges is ignored in
#' directed graphs.
#' @param types An optional vertex type vector to use instead of the
#' \code{type} vertex attribute. You must supply this argument if the graph has
#' no \code{type} vertex attribute.
#' @param attr Either \code{NULL} or a character string giving an edge
#' attribute name. If \code{NULL}, then a traditional incidence matrix is
#' returned. If not \code{NULL} then the values of the given edge attribute are
#' included in the incidence matrix. If the graph has multiple edges, the edge
#' attribute of an arbitrarily chosen edge (for the multiple edges) is
#' included.
#' @param names Logical scalar, if \code{TRUE} and the vertices in the graph
#' are named (i.e. the graph has a vertex attribute called \code{name}), then
#' vertex names will be added to the result as row and column names. Otherwise
#' the ids of the vertices are used as row and column names.
#' @param sparse Logical scalar, if it is \code{TRUE} then a sparse matrix is
#' created, you will need the \code{Matrix} package for this.
#' @return A sparse or dense matrix.
#' @author Gabor Csardi \email{csardi.gabor@@gmail.com}
#' @seealso \code{\link{graph_from_incidence_matrix}} for the opposite operation.
#' @export
#' @keywords graphs
#' @examples
#' 
#' g <- bipartite_graph( c(0,1,0,1,0,0), c(1,2,2,3,3,4) )
#' as_incidence_matrix(g)
#' 
as_incidence_matrix <- function(graph, types=NULL, attr=NULL,
                          names=TRUE, sparse=FALSE) {
  # Argument checks
  if (!is_igraph(graph)) { stop("Not a graph object") }
  if (is.null(types) && "type" %in% vertex_attr_names(graph)) { 
    types <- V(graph)$type 
  } 
  if (!is.null(types)) { 
    types <- as.logical(types) 
  } else { 
    stop("Not a bipartite graph, supply `types' argument") 
  }
  
  names <- as.logical(names)
  sparse <- as.logical(sparse)
  
  if (sparse) {
    get.incidence.sparse(graph, types=types, names=names, attr=attr)
  } else {
    get.incidence.dense(graph, types=types, names=names, attr=attr)
  }
}

#' @rdname graph_from_data_frame
#' @param x An igraph object.
#' @param what Character constant, whether to return info about vertices,
#' edges, or both. The default is \sQuote{edges}.
#' @export

as_data_frame <- function(x, what=c("edges", "vertices", "both")) {

  if (!is_igraph(x)) { stop("Not a graph object") }
  what <- igraph.match.arg(what)

  if (what %in% c("vertices", "both")) {
    ver <- .Call("R_igraph_mybracket2", x, 9L, 3L, PACKAGE="igraph")
    class(ver) <- "data.frame"
    rn <- if (is_named(x)) { V(x)$name } else { seq_len(vcount(x)) }
    rownames(ver) <- rn
  }

  if (what %in% c("edges", "both")) {
    el <- as_edgelist(x)
    edg <- c(list(from=el[,1]), list(to=el[,2]),
             .Call("R_igraph_mybracket2", x, 9L, 4L, PACKAGE="igraph"))
    class(edg) <- "data.frame"
    rownames(edg) <- seq_len(ecount(x))
  }
  
  if (what=="both") {
    list(vertices=ver, edges=edg)
  } else if (what=="vertices") {
    ver
  } else {
    edg
  }
}
