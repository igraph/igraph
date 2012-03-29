
#   IGraph R package
#   Copyright (C) 2005  Gabor Csardi <csardi@rmki.kfki.hu>
#   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
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

V <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  vc <- vcount(graph)
  if (vc == 0) {
    res <- numeric()
  } else {
    res <- 0:(vc-1)
  }
  class(res) <- "igraph.vs"
  ne <- new.env()
  assign("graph", graph, envir=ne)
  attr(res, "env") <- ne
  res
}

E <- function(graph, P=NULL, path=NULL, directed=TRUE) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  if (!is.null(P) && !is.null(path)) {
    stop("Cannot give both `P' and `path' at the same time")
  }
  
  if (is.null(P) && is.null(path)) {  
    ec <- ecount(graph)
    if (ec == 0) {
      res <- numeric()
    } else {
      res <- 0:(ec-1)
    }
  } else if (!is.null(P)) {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
    res <- .Call("R_igraph_es_pairs", graph, as.numeric(P),
                 as.logical(directed),
                 PACKAGE="igraph0")
  } else {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
    res <- .Call("R_igraph_es_path", graph, as.numeric(path),
                 as.logical(directed),
                 PACKAGE="igraph0")
  }
    
  class(res) <- "igraph.es"
  ne <- new.env()
  assign("graph", graph, envir=ne)
  attr(res, "env") <- ne
  res
}

"[.igraph.vs" <- function(x, i) {
  i <- substitute(i)
  if (is.numeric(i) || is.integer(i)) {
    # simple indexing by vertex ids
    res <- i[ i %in% x ]
    attributes(res) <- attributes(x)
  } else if (is.logical(i)) {
    # simple indexing by logical vector
    res <- as.numeric(x) [ i ]
    attributes(res) <- attributes(x)
  } else if (is.character(i)) {
    res <- as.igraph.vs(get("graph", attr(x, "env")), i)
    attributes(res) <- attributes(x)
  } else {
    # language expression, we also do attribute based indexing
    graph <- get("graph", attr(x, "env"))
    nei <- function(v, mode=c("all", "in", "out", "total")) {
      ## TRUE iff the vertex is a neighbor (any type)
      ## of at least one vertex in v
      mode <- igraph.match.arg(mode)
      mode <- switch(mode, "out"=1, "in"=2, "all"=3, "total"=3)

      if (is.logical(v)) {
        v <- which(v)
      }
      on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
      tmp <- .Call("R_igraph_vs_nei", graph, x, as.igraph.vs(graph, v),
                   as.numeric(mode),
                   PACKAGE="igraph0")
      tmp[as.numeric(x)+1]
    }
    innei <- function(v, mode=c("in", "all", "out", "total")) {
      nei(v, mode)
    }
    outnei <- function(v, mode=c("out", "all", "in", "total")) {
      nei(v, mode)
    }
    adj <- function(e) {
      ## TRUE iff the vertex (in the vs) is adjacent
      ## to at least one edge in e
      if (is.logical(e)) {
        e <- which(e)
      }
      on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
      tmp <- .Call("R_igraph_vs_adj", graph, x, as.igraph.es(e), as.numeric(3),
                   PACKAGE="igraph0")
      tmp[as.numeric(x)+1]
    }
    from <- function(e) {
      ## TRUE iff the vertex is the source of at least one edge in e
      if (is.logical(e)) {
        e <- which(e)
      }
      on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
      tmp <- .Call("R_igraph_vs_adj", graph, x, as.igraph.es(e), as.numeric(1),
                   PACKAGE="igraph0")
      tmp[as.numeric(x)+1]
    }
    to <- function(e) {
      ## TRUE iff the vertex is the target of at least one edge in e
      if (is.logical(e)) {
        e <- which(e)
      }
      on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
      tmp <- .Call("R_igraph_vs_adj", graph, x, as.igraph.es(e), as.numeric(2),
                   PACKAGE="igraph0")
      tmp[as.numeric(x)+1]
    }
    i <- eval(i, envir=c(graph[[9]][[3]], nei=nei, innei=innei,
                   outnei=outnei, adj=adj, from=from, to=to),
              enclos=parent.frame())
    if (is.numeric(i) || is.integer(i)) {
      i <- as.numeric(i)
      res <- i[ i %in% x ]
      attributes(res) <- attributes(x)
    } else if (is.logical(i)) {
      res <- as.numeric(x) [ i ]
      attributes(res) <- attributes(x)
    } else if (is.character(i)) {
      res <- as.igraph.vs(get("graph", attr(x, "env")), i)
      attributes(res) <- attributes(x)
    } else {
      stop("invalid indexing of vertex seq")
    }
  }

  res
}

"[.igraph.es" <- function(x, i) {
  i <- substitute(i)
  if (is.numeric(i) || is.integer(i)) {
    # simple indexing by vertex ids
    res <- i[ i %in% x ]
    attributes(res) <- attributes(x)    
  } else if (is.logical(i)) {
    # simple indexing by a logical vector
    res <- as.numeric(x) [ i ]
    attributes(res) <- attributes(x)
  } else {
    # language expression, we also do attribute based indexing
    graph <- get("graph", attr(x, "env"))
    i <- substitute(i)
    adj <- function(v) {
      ## TRUE iff the edge is adjacent to at least one vertex in v
      on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
      tmp <- .Call("R_igraph_es_adj", graph, x, as.igraph.vs(graph, v), as.numeric(3),
                   PACKAGE="igraph0")
      tmp[ as.numeric(x)+1 ]
    }
    from <- function(v) {
      ## TRUE iff the edge originates from at least one vertex in v
      on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
      tmp <- .Call("R_igraph_es_adj", graph, x, as.igraph.vs(graph, v), as.numeric(1),
                   PACKAGE="igraph0")
      tmp[ as.numeric(x)+1 ]      
    }
    to <- function(v) {
      ## TRUE iff the edge points to at least one vertex in v
      on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
      tmp <- .Call("R_igraph_es_adj", graph, x, as.igraph.vs(graph, v), as.numeric(2),
                   PACKAGE="igraph0")
      tmp[ as.numeric(x)+1 ]
    }
    i <- eval(i, envir=c(graph[[9]][[4]],
                   adj=adj, from=from, to=to,
                   .igraph.from=list(graph[[3]][ as.numeric(x)+1 ]),
                   .igraph.to=list(graph[[4]][as.numeric(x)+1]),
                   .igraph.graph=list(graph),
                   `%--%`=`%--%`, `%->%`=`%->%`, `%<-%`=`%<-%`),
              enclos=parent.frame())
    if (is.numeric(i) || is.integer(i)) {
      i <- as.numeric(i)
      res <- i[ i %in% x ]
      attributes(res) <- attributes(x)
    } else if (is.logical(i)) {
      res <- as.numeric(x) [ i ]
      attributes(res) <- attributes(x)
    } else {
      stop("invalid indexing of edge seq")
    }
  }
  
  res
} 

"%--%" <- function(f, t) {
  from <- get(".igraph.from", parent.frame())
  to <- get(".igraph.to", parent.frame())
  (from %in% f & to %in% t) | (to %in% f & from %in% t)
}

"%->%" <- function(f, t) {
  from <- get(".igraph.from", parent.frame())
  to <- get(".igraph.to", parent.frame())
  graph <- get(".igraph.graph", parent.frame())
  if (is.directed(graph)) {
    from %in% f & to %in% t
  } else {
    (from %in% f & to %in% t) | (to %in% f & from %in% t)
  }
}

"%<-%" <- function(t, value) {
  from <- get(".igraph.from", parent.frame())
  to <- get(".igraph.to", parent.frame())
  graph <- get(".igraph.graph", parent.frame())
  if (is.directed(graph)) {
    from %in% value & to %in% t
  } else {
    (from %in% value & to %in% t) | (to %in% value & from %in% t)
  }
}

"[<-.igraph.vs" <- function(x, i, value) {
  if (! "name"  %in% names(attributes(value)) ||
      ! "value" %in% names(attributes(value))) {
    stop("invalid indexing")
  }
  value
}

"[<-.igraph.es" <- function(x, i, value) {
  if (! "name"  %in% names(attributes(value)) ||
      ! "value" %in% names(attributes(value))) {
    stop("invalid indexing")
  }
  value
}  

"$.igraph" <- function(x, name) {
  get.graph.attribute(x, name)
}

"$<-.igraph" <- function(x, name, value) {
  set.graph.attribute(x, name, value)
}
  
"$.igraph.vs" <- function(x, name) {
  get.vertex.attribute(get("graph", attr(x, "env")), name, x)
}

"$.igraph.es" <- function(x, name) {
  get.edge.attribute(get("graph", attr(x, "env")), name, x)
}

"$<-.igraph.vs" <- function(x, name, value) {
  attr(x, "name") <- name
  attr(x, "value") <- value
  x
}

"$<-.igraph.es" <- function(x, name, value) {
  attr(x, "name") <- name
  attr(x, "value") <- value
  x
}

"V<-" <- function(x, value) {
  if (!is.igraph(x)) {
    stop("Not a graph object")
  }
  if (! "name"  %in% names(attributes(value)) ||
      ! "value" %in% names(attributes(value))) {
    stop("invalid indexing")
  }
  set.vertex.attribute(x, attr(value, "name"), index=value,
                       value=attr(value, "value"))
}

"E<-" <- function(x, path=NULL, P=NULL, directed=NULL, value) {
  if (!is.igraph(x)) {
    stop("Not a graph object")
  }
  if (! "name"  %in% names(attributes(value)) ||
      ! "value" %in% names(attributes(value))) {
    stop("invalid indexing")
  }
  set.edge.attribute(x, attr(value, "name"), index=value,
                     value=attr(value, "value"))
}

print.igraph.vs <- function(x, ...) {
  cat("Vertex sequence:\n")
  graph <- get("graph", attr(x, "env"))
  x <- as.numeric(x)
  if ("name" %in% list.vertex.attributes(graph)) {
    x <- V(graph)$name[x+1]
  }
  print(x)
}

print.igraph.es <- function(x, ...) {
  cat("Edge sequence:\n")
  graph <- get("graph", attr(x, "env"))
  if (is.directed(graph)) {
    arrow <- "->"
  } else {
    arrow <- "--"
  }
  x <- as.numeric(x)
  el <- get.edges(graph, x)
  if ("name" %in% list.vertex.attributes(graph)) {
    el <- matrix(V(graph)$name[el+1], ncol=2)
  }
  tab <- data.frame(e=paste(sep="", "[", x, "]"), row.names="e")
  if (is.numeric(el)) { w <- nchar(max(el)) } else { w <- max(nchar(el)) }
  tab[" "] <- paste(format(el[,1], width=w), arrow, format(el[,2], width=w))
  print(tab)
}

# these are internal

as.igraph.vs <- function(graph, v) {
  if (is.character(v) && "name" %in% list.vertex.attributes(graph)) {
    v <- as.numeric(match(v, V(graph)$name)-1)
    if (any(is.na(v))) {
      stop("Invalid vertex names")
    }
    v
  } else {
    as.numeric(v)
  }
}

as.igraph.es <- as.numeric
