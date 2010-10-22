
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

`[.igraph` <- function(x, i, j, ..., directed=TRUE, weighted=TRUE, 
                       multi=FALSE, edges=NULL, simplify=TRUE, unique=TRUE,
                       drop=TRUE) {
  if (!is.null(edges) && (!missing(i) || !missing(j))) {
    stop("Vertices cannot be given together with 'edges'")
  }
  if (!missing(i) && !missing(j)) {
    ## Query about connections
    weighted <- weighted && "weight" %in% list.edge.attributes(x)
    if (length(i) != length(j)) {
      stop("'i' and 'j' must have the same length")
    }
    res <- get.edge.ids(x, rbind(i,j), directed=directed, error=FALSE,
                        multi=multi)
    if (weighted) {
      res <- rep(NA, length(i))
      res[res!=0] <- get.edge.attribute(x, "weight", res[res!=0])
    } else {
      res <- as.logical(res)
    }
  } else if (!is.null(edges)) {
    ## Query incident vertices for a set of edges
    on.exit(.Call("R_igraph_finalizer", PACKAGE="igraph"))
    res <- .Call("R_igraph_edges", x, as.igraph.es(x, edges)-1,
                 PACKAGE="igraph")+1
    if (!simplify) {
      ## TODO: rewrite this in C to make it faster
      res <- split(res, rep(seq_along(edges), each=2))
    } else if (unique) {
      res <- unique(res)
    }
  } else {
    ## Query adjacenct vertices for some vertices
    ## TODO: do this with one C call
    if (missing(i)) {
      mode <- if (directed) "in" else "all"
      v <- j
    } else {
      mode <- if (directed) "out" else "all"
      v <- i
    }
    res <- lapply(v, neighbors, graph=x, mode=mode)
    if (simplify) {
      res <- unlist(res)
      if (unique) { res <- unique(res) }
    }
  }
  res
}

`[[.igraph` <- function(x, i, j, ..., directed=TRUE, multi=FALSE,
                        simplify=TRUE, unique=TRUE, exact=TRUE) {

  if (!missing(i) && !missing(j)) {
    ## Query edge ids connection some vertices
    if (length(i) != length(j)) {
      stop("'i' and 'j' must have the same length")
    }
    res <- get.edge.ids(x, rbind(i,j), directed=directed, error=FALSE,
                        multi=multi)
    res[res==0] <- NA
  } else {
    ## Incident edges for some vertices
    ## TODO: do this with one C call
    if (missing(i)) {
      mode <- if (directed) "in" else "all"
      v <- j
    } else {
      mode <- if (directed) "out" else "all"
      v <- i
    }
    res <- lapply(v, incident, graph=x, mode=mode)
    if (simplify) {
      res <- unlist(res)
      if (unique) { res <- unique(res) }
    }
  }
  res
}

`[<-.igraph` <- function(x, i, j, ..., add=FALSE, value) {

  if (!is.null(value) && !is.numeric(value) && !is.logical(value)) {
    stop("The new value should be NULL, numeric or logical")
  }
  if (is.logical(value) && length(value) != 1) {
    stop("If new value if logical, then it should be of length 1")
  }
  if (is.logical(value) && is.na(value)) {
    stop("Logical value cannot be NA")
  }
  
  if (is.null(value) || (is.logical(value) && !value)) {
    ## Deletion, get the edges to be deleted
    todel <- x[[i, j, ...]]
    x <- delete.edges(x, todel)
  } else {
    ## Not deletion, can be addition or update or mixed
    add <- add | ! x[i,j,...]
    if (any(add)) {
      if (is.logical(value)) {
        ## Not weighted
        toadd <- as.igraph.vs(x, as.vector(rbind(i[add],j[add])))
        x <- add.edges(x, toadd)
      } else {
        ## Weighted
        if (length(value) != 1 && length(value) != length(i)) {
          stop("Invalid 'value' length, should be 1 or the same as the ",
               "number of edges")
        }
        if (! "weight" %in% list.edge.attributes(x)) {
          x <- set.edge.attribute(x, "weight", value=NA)
        }
        toadd <- as.igraph.vs(x, as.vector(rbind(i[add],j[add])))
        x <- add.edges(x, toadd, weight=value)
      }
    }
    if (any(!add) && is.numeric(value)) {
      ## Update weights
      toupd <- x[[i[!add],j[!add],...]]
      x <- set.edge.attribute(x, "weight", toupd, value[!add])
    }
  }
  x
}
