
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
# Constructors
###################################################################

update_es_ref <- update_vs_ref <- function(graph) {
  env <- get_vs_ref(graph)
  if (!is.environment(env)) {
    stop("This graph was created by an old(er) igraph version.\n",
         "  Call upgrade_graph() on it to use with the current igraph version")
  }

  assign("me", graph, envir = env)
}

get_es_ref <- get_vs_ref <- function(graph) {
  .Call("R_igraph_mybracket", graph, 10L, PACKAGE = "igraph")
}

get_es_graph <- get_vs_graph <- function(seq) {
  weak_ref_key(attr(seq, "env"))$me
}

has_es_graph <- has_vs_graph <- function(seq) {
  !is.null(weak_ref_key(attr(seq, "env")))
}

#' Decide if two graphs are identical
#'
#' This is similar to \code{identical} in the \code{base} package,
#' but ignores the mutable piece of igraph objects, that might be
#' different, even if the two graphs are identical.
#'
#' @param g1,g2 The two graphs
#' @return Logical scalar
#' @export

identical_graphs <- function(g1, g2) {
  stopifnot(is_igraph(g1), is_igraph(g2))
  .Call("R_igraph_identical_graphs", g1, g2, PACKAGE = "igraph");
}

#' @export

V <- function(graph) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  update_vs_ref(graph)

  res <- seq_len(vcount(graph))
  if (is_named(graph)) names(res) <- vertex_attr(graph)$name
  class(res) <- "igraph.vs"
  attr(res, "env") <- make_weak_ref(get_vs_ref(graph), NULL)
  res
}

create_vs <- function(graph, idx, na_ok = FALSE) {
  if (na_ok) idx <- ifelse(idx < 1 | idx > gorder(graph), NA, idx)
  V(graph)[idx]
}

#' @export

E <- function(graph, P=NULL, path=NULL, directed=TRUE) {
  if (!is_igraph(graph)) {
    stop("Not a graph object")
  }

  update_es_ref(graph)

  if (!is.null(P) && !is.null(path)) {
    stop("Cannot give both `P' and `path' at the same time")
  }
  
  if (is.null(P) && is.null(path)) {  
    ec <- ecount(graph)
    res <- seq_len(ec)
  } else if (!is.null(P)) {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    res <- .Call("R_igraph_es_pairs", graph, as.igraph.vs(graph, P)-1,
                 as.logical(directed),
                 PACKAGE="igraph")+1
  } else {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    res <- .Call("R_igraph_es_path", graph, as.igraph.vs(graph, path)-1,
                 as.logical(directed),
                 PACKAGE="igraph")+1
  }

  if ("name" %in% edge_attr_names(graph)) {
    names(res) <- edge_attr(graph)$name[res]
  } else if (is_named(graph)) {
    el <- as_edgelist(graph)
    names(res) <- paste(el[,1], el[,2], sep = "|")
  }
  
  class(res) <- "igraph.es"
  attr(res, "env") <- make_weak_ref(get_es_ref(graph), NULL)
  res
}

create_es <- function(graph, idx, na_ok = FALSE) {
  if (na_ok) idx <- ifelse(idx < 1 | idx > gsize(graph), NA, idx)
  E(graph)[idx]
}

## TODO: remove this? What is the point?
#' @method "[[" igraph.vs
#' @export

"[[.igraph.vs" <- function(x, i) {
  if (length(i) != 1) {
    stop("Invalid `[[` indexing, need single vertex")
  }
  res <- x[i]
  attr(res, "single") <- TRUE
  res
}

simple_vs_index <- function(x, i) {
  res <- unclass(x)[i]
  if (anyNA(res)) stop('Unknown vertex selected')
  res
}

#' @method "[" igraph.vs
#' @export

"[.igraph.vs" <- function(x, i) {
  graph <- get_vs_graph(x)
  if (is.null(graph)) {
    res <- simple_vs_index(x, i)
  } else {
    i <- substitute(i)
    nei <- function(v, mode=c("all", "in", "out", "total")) {
      ## TRUE iff the vertex is a neighbor (any type)
      ## of at least one vertex in v
      mode <- igraph.match.arg(mode)
      mode <- switch(mode, "out"=1, "in"=2, "all"=3, "total"=3)

      if (is.logical(v)) {
        v <- which(v)
      }
      on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
      tmp <- .Call("R_igraph_vs_nei", graph, x, as.igraph.vs(graph, v)-1,
                   as.numeric(mode),
                   PACKAGE="igraph")
      tmp[as.numeric(x)]
    }
    innei <- function(v, mode=c("in", "all", "out", "total")) {
      nei(v, mode)
    }
    outnei <- function(v, mode=c("out", "all", "in", "total")) {
      nei(v, mode)
    }
    inc <- adj <- function(e) {
      ## TRUE iff the vertex (in the vs) is incident
      ## to at least one edge in e
      if (is.logical(e)) {
        e <- which(e)
      }
      on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
      tmp <- .Call("R_igraph_vs_adj", graph, x, as.igraph.es(graph, e)-1,
                   as.numeric(3),
                   PACKAGE="igraph")
      tmp[as.numeric(x)]
    }
    from <- function(e) {
      ## TRUE iff the vertex is the source of at least one edge in e
      if (is.logical(e)) {
        e <- which(e)
      }
      on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
      tmp <- .Call("R_igraph_vs_adj", graph, x, as.igraph.es(graph, e)-1,
                   as.numeric(1),
                   PACKAGE="igraph")
      tmp[as.numeric(x)]
    }
    to <- function(e) {
      ## TRUE iff the vertex is the target of at least one edge in e
      if (is.logical(e)) {
        e <- which(e)
      }
      on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
      tmp <- .Call("R_igraph_vs_adj", graph, x, as.igraph.es(graph, e)-1,
                   as.numeric(2),
                   PACKAGE="igraph")
      tmp[as.numeric(x)]
    }
    i <- eval(i, envir=c(.Call("R_igraph_mybracket2", graph, 9L, 3L,
                   PACKAGE="igraph"), nei=nei, innei=innei,
                   outnei=outnei, adj=adj, inc=inc, from=from, to=to),
              enclos=parent.frame())
    res <- simple_vs_index(x, i)
  }
  attr(res, "env") <- attr(x, "env")
  class(res) <- class(x)
  res
}

## TODO: remove this? What is it for?
#' @method "[[" igraph.es
#' @export

"[[.igraph.es" <- function(x, i) {
  if (length(i) != 1) {
    stop("Invalid `[[` indexing, need single edge")
  }
  res <- x[i]
  attr(res, "single") <- TRUE
  res
}

simple_es_index <- function(x, i) {
  res <- unclass(x)[i]
  if (anyNA(res)) stop('Unknown edge selected')
  res
}

#' @method "[" igraph.es
#' @export

"[.igraph.es" <- function(x, i) {
  i <- substitute(i)
  graph <- get_es_graph(x)
  if (is.null(graph)) {
    res <- simple_es_index(x, i)
  } else {
    ## language expression, we also do attribute based indexing
    i <- substitute(i)
    inc <- adj <- function(v) {
      ## TRUE iff the edge is incident to at least one vertex in v
      on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
      tmp <- .Call("R_igraph_es_adj", graph, x, as.igraph.vs(graph, v)-1,
                   as.numeric(3),
                   PACKAGE="igraph")
      tmp[ as.numeric(x) ]
    }
    from <- function(v) {
      ## TRUE iff the edge originates from at least one vertex in v
      on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
      tmp <- .Call("R_igraph_es_adj", graph, x, as.igraph.vs(graph, v)-1,
                   as.numeric(1),
                   PACKAGE="igraph")
      tmp[ as.numeric(x) ]      
    }
    to <- function(v) {
      ## TRUE iff the edge points to at least one vertex in v
      on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
      tmp <- .Call("R_igraph_es_adj", graph, x, as.igraph.vs(graph, v)-1,
                   as.numeric(2),
                   PACKAGE="igraph")
      tmp[ as.numeric(x) ]
    }
    i <- eval(i, envir=c(.Call("R_igraph_mybracket2", graph, 9L, 4L,
                   PACKAGE="igraph"),
                   inc=inc, adj=adj, from=from, to=to,
                   .igraph.from=list(.Call("R_igraph_mybracket",
                     graph, 3L, PACKAGE="igraph")[ as.numeric(x) ]),
                   .igraph.to=list(.Call("R_igraph_mybracket",
                     graph, 4L, PACKAGE="igraph")[as.numeric(x)]),
                   .igraph.graph=list(graph),
                   `%--%`=`%--%`, `%->%`=`%->%`, `%<-%`=`%<-%`),
              enclos=parent.frame())
    res <- simple_es_index(x, i)
  }
  attr(res, "env") <- attr(x, "env")
  class(res) <- class(x)
  res
} 

#' @export

"%--%" <- function(f, t) {  
  from <- get(".igraph.from", parent.frame())
  to <- get(".igraph.to", parent.frame())
  graph <- get(".igraph.graph", parent.frame())
  f <- as.igraph.vs(graph, f)-1
  t <- as.igraph.vs(graph, t)-1
  (from %in% f & to %in% t) | (to %in% f & from %in% t)
}

#' @export

"%->%" <- function(f, t) {
  from <- get(".igraph.from", parent.frame())
  to <- get(".igraph.to", parent.frame())
  graph <- get(".igraph.graph", parent.frame())
  f <- as.igraph.vs(graph, f)-1
  t <- as.igraph.vs(graph, t)-1
  if (is_directed(graph)) {
    from %in% f & to %in% t
  } else {
    (from %in% f & to %in% t) | (to %in% f & from %in% t)
  }
}

#' @export

"%<-%" <- function(t, value) {
  from <- get(".igraph.from", parent.frame())
  to <- get(".igraph.to", parent.frame())
  graph <- get(".igraph.graph", parent.frame())
  value <- as.igraph.vs(graph, value)-1
  t <- as.igraph.vs(graph, t)-1
  if (is_directed(graph)) {
    from %in% value & to %in% t
  } else {
    (from %in% value & to %in% t) | (to %in% value & from %in% t)
  }
}

#' @method "[[<-" igraph.vs
#' @export

"[[<-.igraph.vs" <- function(x, i, value) {
  if (! "name"  %in% names(attributes(value)) ||
      ! "value" %in% names(attributes(value))) {
    stop("invalid indexing")
  }
  value
}

#' @method "[<-" igraph.vs
#' @export

"[<-.igraph.vs" <-  `[[<-.igraph.vs`

#' @method "[[<-" igraph.es
#' @export

"[[<-.igraph.es" <- function(x, i, value) {
  if (! "name"  %in% names(attributes(value)) ||
      ! "value" %in% names(attributes(value))) {
    stop("invalid indexing")
  }
  value
}  

#' @method "[<-" igraph.es
#' @export

"[<-.igraph.es" <-  `[[<-.igraph.es`

#' @method "$" igraph
#' @export

"$.igraph" <- function(x, name) {
  graph_attr(x, name)
}

#' @method "$<-" igraph
#' @export

"$<-.igraph" <- function(x, name, value) {
  set_graph_attr(x, name, value)
}

#' @method "$" igraph.vs
#' @export

"$.igraph.vs" <- function(x, name) {
  graph <- get_vs_graph(x)
  if (is.null(graph)) stop("Graph is unknown")
  res <- vertex_attr(graph, name, x)
  if ("single" %in% names(attributes(x)) && attr(x, "single")) {
    res[[1]]
  } else {
    res
  }
}
#' @method "$" igraph.es
#' @export

"$.igraph.es" <- function(x, name) {
  graph <- get_es_graph(x)
  if (is.null(graph)) stop("Graph is unknown")
  res <- edge_attr(graph, name, x)
  if ("single" %in% names(attributes(x)) && attr(x, "single")) {
    res[[1]]
  } else {
    res
  }
}

#' @method "$<-" igraph.vs
#' @export

"$<-.igraph.vs" <- function(x, name, value) {
  attr(x, "name") <- name
  attr(x, "value") <- value
  x
}

#' @method "$<-" igraph.es
#' @export

"$<-.igraph.es" <- function(x, name, value) {
  attr(x, "name") <- name
  attr(x, "value") <- value
  x
}

#' @export

"V<-" <- function(x, value) {
  if (!is_igraph(x)) {
    stop("Not a graph object")
  }
  if (! "name"  %in% names(attributes(value)) ||
      ! "value" %in% names(attributes(value))) {
    stop("invalid indexing")
  }
  i_set_vertex_attr(x, attr(value, "name"), index=value,
                    value=attr(value, "value"), check = FALSE)
}

#' @export

"E<-" <- function(x, path=NULL, P=NULL, directed=NULL, value) {
  if (!is_igraph(x)) {
    stop("Not a graph object")
  }
  if (! "name"  %in% names(attributes(value)) ||
      ! "value" %in% names(attributes(value))) {
    stop("invalid indexing")
  }
  i_set_edge_attr(x, attr(value, "name"), index=value,
                  value=attr(value, "value"), check = FALSE)
}

#' @method print igraph.vs
#' @export

print.igraph.vs <- function(x, ...) {
  graph <- get_vs_graph(x)
  if (is.null(graph)) stop("Graph is unknown")
  cat("+ vertices", if (is_named(graph)) ", named" else "", ":\n", sep = "")
  n <- names(x)
  x <- as.numeric(x)
  if ("name" %in% vertex_attr_names(graph)) {
    x <- V(graph)$name[x]
  }
  names(x) <- n
  print(x, quote = FALSE)
  invisible(x)
}

#' @method print igraph.es
#' @export

print.igraph.es <- function(x, ...) {
  graph <- get_es_graph(x)
  if (is.null(graph)) stop("Graph is unknown")
  .print.edges.compressed(graph, x, names = TRUE)
  invisible(x)
}

# these are internal

as.igraph.vs <- function(graph, v, na.ok=FALSE) {
  if (inherits(v, "igraph.vs") && has_vs_graph(v)) {
    if (address(graph) != address(get_vs_graph(v))) {
      stop("Cannot use a vertex sequence from another graph.")
    }
  }
  if (is.character(v) && "name" %in% vertex_attr_names(graph)) {
    v <- as.numeric(match(v, V(graph)$name))
    if (!na.ok && any(is.na(v))) {
      stop("Invalid vertex names")
    }
    v
  } else {
    if (is.logical(v)) {
      res <- as.vector(V(graph))[v]
    } else if (is.numeric(v) && any(v<0)){
      res <- as.vector(V(graph))[v]
    } else {
      res <- as.numeric(v)
    }
    if (!na.ok && any(is.na(res))) {
      stop("Invalid vertex name(s)")
    }
    res
  }
}

as.igraph.es <- function(graph, e) {
  if (inherits(e, "igraph.es") && has_es_graph(e)) {
    if (address(graph) != address(get_es_graph(e))) {
      stop("Cannot use an edge sequence from another graph.")
    }
  }
  if (is.character(e)) {
    Pairs <- grep("|", e, fixed=TRUE)
    Names <- if (length(Pairs)==0) seq_along(e) else -Pairs
    res <- numeric(length(e))

    ## Based on vertex ids/names
    if (length(Pairs)!=0) {
      vv <- strsplit(e[Pairs], "|", fixed=TRUE)
      vl <- sapply(vv, length)
      if (any(vl != 2)) {
        stop("Invalid edge name: ", e[Pairs][vl!=2][1])
      }
      vp <- unlist(vv)
      if (! "name" %in% vertex_attr_names(graph)) {
        vp <- as.numeric(vp)
      }
      res[Pairs] <- get.edge.ids(graph, vp)
    }

    ## Based on edge ids/names
    if (length(Names) != 0) {
      if ("name" %in% edge_attr_names(graph)) {
        res[Names] <- as.numeric(match(e[Names], E(graph)$name))
      } else {
        res[Names] <- as.numeric(e[Names])
      }
    }
    
  } else {
    res <- as.numeric(e)
  }
  if (any(is.na(res))) {
    stop("Invalid edge names")
  }
  res
}
