
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
  res <- seq_len(vc)
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
    i <- eval(i, envir=c(unclass(graph)[[9]][[3]], nei=nei, innei=innei,
                   outnei=outnei, adj=adj, inc=inc, from=from, to=to),
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
  } else if (is.character(i)) {
    res <- as.igraph.es(get("graph", attr(x, "env")), i)
    attributes(res) <- attributes(x)
  } else {
    # language expression, we also do attribute based indexing
    graph <- get("graph", attr(x, "env"))
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
    i <- eval(i, envir=c(unclass(graph)[[9]][[4]],
                   inc=inc, adj=adj, from=from, to=to,
                   .igraph.from=list(unclass(graph)[[3]][ as.numeric(x) ]),
                   .igraph.to=list(unclass(graph)[[4]][as.numeric(x)]),
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
  graph <- get(".igraph.graph", parent.frame())
  f <- as.igraph.vs(graph, f)-1
  t <- as.igraph.vs(graph, t)-1
  (from %in% f & to %in% t) | (to %in% f & from %in% t)
}

"%->%" <- function(f, t) {
  from <- get(".igraph.from", parent.frame())
  to <- get(".igraph.to", parent.frame())
  graph <- get(".igraph.graph", parent.frame())
  f <- as.igraph.vs(graph, f)-1
  t <- as.igraph.vs(graph, t)-1
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
  value <- as.igraph.vs(graph, value)-1
  t <- as.igraph.vs(graph, t)-1
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
    x <- V(graph)$name[x]
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
    el <- matrix(V(graph)$name[el], nc=2)
  }
  tab <- data.frame(e=paste(sep="", "[", x, "]"), row.names="e")
  if (is.numeric(el)) { w <- nchar(max(el)) } else { w <- max(nchar(el)) }
  tab[" "] <- paste(format(el[,1], width=w), arrow, format(el[,2], width=w))
  print(tab)
}

# these are internal

as.igraph.vs <- function(graph, v) {
  if (is.character(v) && "name" %in% list.vertex.attributes(graph)) {
    v <- as.numeric(match(v, V(graph)$name))
    if (any(is.na(v))) {
      stop("Invalid vertex names")
    }
    v
  } else if (is.logical(v)) {
    as.vector(V(graph))[v]
  } else if (is.numeric(v) && any(v<=0)){
    as.vector(V(graph))[v]
  } else {
    as.numeric(v)
  }
}

as.igraph.es <- function(graph, e) {
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
      if (! "name" %in% list.vertex.attributes(graph)) {
        vp <- as.numeric(vp)
      }
      res[Pairs] <- get.edge.ids(graph, vp)
    }

    ## Based on edge ids/names
    if (length(Names) != 0) {
      if ("name" %in% list.edge.attributes(graph)) {
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
