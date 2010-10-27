
## IGraph library.
## Copyright (C) 2010  Gabor Csardi <csardi.gabor@gmail.com>
## Rue de l'Industrie 5, Lausanne 1005, Switzerland

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
## 02110-1301 USA

# Indexing of igraph graphs.
#
# Goals:
# 1. flexible graph manipulation
# 2. to be as close to the usual matrix and adjacency list semantics,
#    as possible
# 3. simple
# 4. fast
# 5. orthogonal
#
# Rules:
# - [ is about the existence of the edges.
# - [ can be used for weights as well, if the graph is weighted.
# - [[ is about adjacent vertices, and essentially works as an
#   adjacency list.
#
# Use cases:
# - G[1,2]      is there an edge from vertex 1 to vertex 2?
# - G[1,1:3]    are there edges from vertex 1 to vertices 1:3?
# - G[1:2,1:3]  are there adges from vertices 1:2 to vertices 1:3?
#               this returns a (possibly sparse) matrix.
# - G[degree(G)==0,1:4]
#               logical vectors work
# - G[1,-1]     negative indices work
#
# - G[[1,]]     adjacent vertices of 1
# - G[[,1]]     adjacent predessors of 1
# - G[[degree(G),]]
#               logical vectors work
# - G[[-1,]]    negative indices work
#
# - G[1,2,attr="value"]
#               query an edge attribute
# - G[1:3,2,eid=TRUE]
#               create an edge sequence


`[.igraph` <- function(x, i, j, ..., from, to,
                       sparse=getIgraphOpt("sparsematrices"),
                       edges=FALSE, drop=TRUE,
                       attr=if (is.weighted(x)) "weight" else NULL) {
  ## TODO: make it faster, don't need the whole matrix usually

  ################################################################
  ## Argument checks
  if ((!missing(from) || !missing(to)) &&
      (!missing(i)    || !missing(j))) {
    stop("Cannot give 'from'/'to' together with regular indices")
  }
  if ((!missing(from) &&  missing(to)) ||
      ( missing(from) && !missing(to))) {
    stop("Cannot give 'from'/'to' without the other")
  }
  if (!missing(from)) {
    if ((!is.numeric(from) && !is.character(from)) || any(is.na(from))) {
      stop("'from' must be a numeric or character vector without NAs")
    }
    if ((!is.numeric(to) && !is.character(to)) || any(is.na(to))) {
      stop("'to' must be a numeric or character vector without NAs")
    }
    if (length(from) != length(to)) {
      stop("'from' and 'to' must have the same length")
    }
  }

  ##################################################################
  
  if (!missing(from)) {
    res <- get.edge.ids(x, rbind(from, to), error=FALSE)
    if (edges) {
      ## nop
    } else if (!is.null(attr)) {
      if (any(res!=0)) {
        res[res!=0] <- get.edge.attribute(x, attr, res[res!=0])
      }
    } else {
      res <- as.logical(res)+0
    }
    res
  } else if (missing(i) && missing(j)) {
    get.adjacency(x, sparse=sparse, attr=attr, eids=edges)
  } else if (missing(j)) {
    get.adjacency(x, sparse=sparse, attr=attr, edges=edges)[j,,drop=drop]
  } else if (missing(i)) {
    get.adjacency(x, sparse=sparse, attr=attr, edges=edges)[,i,drop=drop]
  } else {
    get.adjacency(x, sparse=sparse, attr=attr, edges=edges)[i,j,drop=drop]
  }
}

`[[.igraph` <- function(x, i, j, ..., directed=TRUE,
                        edges=FALSE, exact=TRUE) {
  ## TODO: make it faster, don't need the whole list usually
  getfun <- if (edges) get.adjedgelist else get.adjlist
  if (missing(i) && missing(j)) {
    mode <- if (directed) "out" else "all"
    getfun(x, mode=mode)
  } else if (missing(j)) {
    mode <- if (directed) "out" else "all"
    getfun(x, mode=mode)[i]
  } else if (missing(i)) {
    mode <- if (directed) "in" else "all"
    getfun(x, mode=mode)[j]
  } else {
    mode <- if (directed) "out" else "all"
    i <- as.igraph.vs(x, i)
    j <- as.igraph.vs(x, j)
    if (!edges) {
      lapply(getfun(x, mode=mode)[i], intersect, j)
    } else {
      ee <- get.adjedgelist(x, mode=mode)[i]
      lapply(seq_along(i), function(yy) {
        from <- i[yy]
        el <- get.edges(x, ee[[yy]])
        other <- ifelse(el[,1]==from, el[,2], el[,1])
        ee[[yy]][other %in% j]
      })
      
    }
  }
}

`[<-.igraph` <- function(x, i, j, ..., from, to,
                         attr=if (is.weighted(x)) "weight" else NULL,
                         value) {
  ## TODO: rewrite this in C to make it faster

  ################################################################
  ## Argument checks
  if ((!missing(from) || !missing(to)) &&
      (!missing(i)    || !missing(j))) {
    stop("Cannot give 'from'/'to' together with regular indices")
  }
  if ((!missing(from) &&  missing(to)) ||
      ( missing(from) && !missing(to))) {
    stop("Cannot give 'from'/'to' without the other")
  }
  if (is.null(attr) &&
      (!is.null(value) && !is.numeric(value) && !is.logical(value))) {
    stop("New value should be NULL, numeric or logical")
  }  
  if (is.null(attr) && !is.null(value) && length(value) != 1) {
    stop("Logical or numeric value must be of length 1")
  }
  if (!missing(from)) {
    if ((!is.numeric(from) && !is.character(from)) || any(is.na(from))) {
      stop("'from' must be a numeric or character vector without NAs")
    }
    if ((!is.numeric(to) && !is.character(to)) || any(is.na(to))) {
      stop("'to' must be a numeric or character vector without NAs")
    }
    if (length(from) != length(to)) {
      stop("'from' and 'to' must have the same length")
    }
  }

  ##################################################################

  if (!missing(from)) {    
    if (is.null(value) ||
        (is.logical(value) && !value) ||
        (is.numeric(value) && value==0)) {
      ## Delete edges
      todel <- x[from=from, to=to, ..., edges=TRUE]
      x <- delete.edges(x, todel)
    } else {
      ## Addition or update of an attribute (or both)
      ids <- x[from=from, to=to, ..., edges=TRUE]
      if (any(ids==0)) {
        x <- add.edges(x, rbind(from[ids==0], to[ids==0]))
      }
      if (!is.null(attr)) {
        ids <- x[from=from, to=to, ..., edges=TRUE]
        x <- set.edge.attribute(x, attr, ids, value=value)
      }
    }
  } else if (is.null(value) ||
      (is.logical(value) && !value) ||
      (is.numeric(value) && value==0)) {
    ## Delete edges
    if (missing(i) && missing(j)) {
      todel <- unlist(x[[ ,  , ..., edges=TRUE]])
    } else if (missing(j)) {
      todel <- unlist(x[[i,  , ..., edges=TRUE]])
    } else if (missing(i)) {
      todel <- unlist(x[[ , j, ..., edges=TRUE]])
    } else {
      todel <- unlist(x[[i, j, ..., edges=TRUE]])
    }
    x <- delete.edges(x, todel)
  } else {
    ## Addition or update of an attribute (or both)
    i <- if (missing(i)) as.numeric(V(x)) else as.igraph.vs(x, i)
    j <- if (missing(j)) as.numeric(V(x)) else as.igraph.vs(x, j)
    if (length(i) != 0 && length(j) != 0) {
      ## Existing edges, and their endpoints
      exe <- x[[i, j, ..., edges=TRUE]]
      exv <- x[[i, j, ...]]
      toadd <- unlist(lapply(seq_along(exv), function(idx) {
        to <- setdiff(j, exv[[idx]])
        if (length(to!=0)) {
          rbind(i[idx], setdiff(j, exv[[idx]]))
        } else {
          numeric()
        }
      }))
      ## Do the changes
      if (is.null(attr)) {
        x <- add.edges(x, toadd)
      } else {
        x <- add.edges(x, toadd, attr=structure(list(value), names=attr))
        toupdate <- unlist(x[[i, j, ..., edges=TRUE]])
        x <- set.edge.attribute(x, attr, toupdate, value)
      }
    }    
  }
  x
}
