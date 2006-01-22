
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
#   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
###################################################################

###################################################################
# Implementation issues.
# There are two types of iterators (this is only true for the R
# interface of course). The first type is a sequence of type
# seq(from, to, by=1). The second type is a simple vector (ie. vector
# of double constants, numeric()).
#
# Both types are represented as a list. The first element of the list
# is always a numeric(3) vector with the following elements:
#   it[[1]][1]  - the actual position of the iterator
#   it[[1]][2]  - the reset (first) position of the iterator
#   it[[1]][3]  - the last position of the iterator
# This common iterator-type independent structure makes it possible to
# implement the most common operations (get, end, next) very
# efficiently, it is not even required to check the type of the
# iterator. See the implementation below.
#
# Sequences are quite simple, (from, to) = it[[1]][2:3].
#
# Vectors are simple as well, they have the whole vector in it[[2]],
# it[[1]][2] is always 1 and it[[1]][3] is the length of the vector.
#
# This is (yet) the only part of the R interface which has some
# built-in knowledge about the structure of the igraph_t object
# (not counting the conversion code in rinterface.c of course).
# It is not good to have such code. But i decided to have because
# (1) it is simple and brief, can be modified easily in case of need
# and (2) otherwise it would be painfully slow. Beacuse of (2) i had
# two choices basically, either leave out iterators from the R
# interface completely or do it in an igraph_t implementation-specific
# way. I considered the latter to be a better choice.

###################################################################
# Constructors
###################################################################

igraph.vs.all <- function(graph) {
  it <- list( c(0,0,vcount(graph)-1,) )
  class(it) <- "igraphvsseq"
  it
}

igraph.es.all <- function(graph) {
  it <- list( c(0,0,ecount(graph)-1,) )
  class(it) <- "igraphesseq"
  it
}

igraph.es.fromorder <- function(graph) {
  .Call("R_igraph_es_fromorder", graph, PACKAGE="igraph")
}

igraph.es.adj <- function(graph, vid, mode="all") {
  if (is.character(mode)) {
    mode <- switch(mode, "out"=1, "in"=2, "all"=3, "total"=3)
  }
  .Call("R_igraph_es_adj", graph, as.numeric(vid), as.numeric(mode),
        PACKAGE="igraph")
}

igraph.vs.adj <- function(graph, vid, mode="all") {
  if (is.character(mode)) {
    mode <- switch(mode, "out"=1, "in"=2, "all"=3, "total"=3)
  }
  .Call("R_igraph_vs_adj", graph, as.numeric(vid), as.numeric(mode),
        PACKAGE="igraph")
}

igraph.vs.vector <- function(graph, v) {
  it <- list( c(1,1,length(v),), as.numeric(v) )
  class(it) <- "igraphvsvector"
  it
}

igraph.es.vector <- function(graph, v) {
  it <- list( c(1,1,length(v),), as.numeric(v) )
  class(it) <- "igraphesvector"
  it
}

###################################################################
# Generic operations
###################################################################

igraph.vs.next <- function(graph, iterator) {
  iterator[[1]][1] <- iterator[[1]][1] + 1
  iterator
}

igraph.vs.end <- function(graph, iterator) {
  iterator[[1]][1] > iterator[[1]][3]
}

# This might be replaced later if it turns out to be too slow
# The same applies to igraph.es.get, igraph.es.from & igraph.es.to

igraph.vs.get <- function(graph, iterator)
  UseMethod("igraph.vs.get", iterator)

igraph.vs.get.igraphvsseq <- function(graph, iterator) {
  iterator[[1]][1]
}

igraph.vs.get.igraphvsvector <- function(graph, iterator) {
  iterator[[2]] [ iterator[[1]][1] ]
}

igraph.vs.reset <- function(graph, iterator) {
  iterator[[1]][1] <- iterator[[1]][2]
}

as.vector.igraphvsseq <- function(x, mode="any") {
  if (x[[1]][2] <= x[[1]][3]) {
    (x[[1]][2]):(x[[1]][3])
  } else {
    numeric()
  }
}

as.vector.igraphvsvector <- function(x, mode="any") {
  x[[2]]
}

igraph.es.next <- function(graph, iterator) {
  iterator[[1]][1] <- iterator[[1]][1] + 1
  iterator
}

igraph.es.end <- function(graph, iterator) {
  iterator[[1]][1] > iterator[[1]][3]
}

igraph.es.get <- function(graph, iterator) 
  UseMethod("igraph.es.get", iterator)

igraph.es.get.igraphesseq <- function(graph, iterator) {
  iterator[[1]][1]
}

igraph.es.get.igraphesvector <- function(graph, iterator) {
  iterator[[2]] [ iterator[[1]][1] ]
}

igraph.es.reset <- function(graph, iterator) {
  iterator[[1]][1] <- iterator[[1]][2]
}

igraph.es.from <- function(graph, iterator)
  UseMethod("igraph.es.from", iterator)

igraph.es.from.igraphesseq <- function(graph, iterator) {
  graph[[3]] [ iterator[[1]][1]+1 ]
}

igraph.es.from.igraphesvector <- function(graph, iterator) {
  graph[[3]] [ iterator[[2]] [ iterator[[1]][1] ]+1 ]
}

igraph.es.to <- function(graph, iterator)
  UseMethod("igraph.es.to", iterator)

igraph.es.to.igraphesseq <- function(graph, iterator) {
  graph[[4]] [ iterator[[1]][1]+1 ]
}

igraph.es.to.igraphesvector <- function(graph, iterator) {
  graph[[4]] [ iterator[[2]] [ iterator[[1]][1] ]+1 ]
}

as.vector.igraphesseq <- function(x, mode="any") {
  if (x[[1]][2] <= x[[1]][3]) {
    (x[[1]][2]):(x[[1]][3])
  } else {
    numeric()
  }
}

as.vector.igraphesvector <- function(x, mode="any") {
  x[[2]]
}

###################################################################
# Iterator shorthands, the should make iterators much more
# comfortable. Look at something like this:
#
# i <- igraph.es.adj(graph, from=c(1,2,3), to=igraph_vs_all(graph))
# while (! i$end) {
#   print(paste(i$e, ":", i$from(g), "->", i$to(g)))
#   i <- i$step
# }
# 
# Isn't it nice?

"$.igraphvsseq" <- function(it, cmd) {
  if (cmd=="end") {
    igraph.vs.end(NULL, it)
  } else if (cmd=="v") {
    igraph.vs.get.igraphvsseq(NULL, it)
  } else if (cmd=="step") {
    igraph.vs.next(NULL, it)
  } else {
    NULL
  }
}

"$.igraphvsvector" <- function(it, cmd) {
  if (cmd=="end") {
    igraph.vs.end(NULL, it)
  } else if (cmd=="v") {
    igraph.vs.get.igraphvsvector(NULL, it)
  } else if (cmd=="step") {
    igraph.vs.next(NULL, it)
  } else {
    NULL
  }
}

"$.igraphesseq" <- function(it, cmd) {
  if (cmd=="end") {
    igraph.es.end(NULL, it)
  } else if (cmd=="e") {
    igraph.es.get.igraphesseq(NULL, it)
  } else if (cmd=="from") {
    function(g) { igraph.es.from.igraphesseq(g, it) }
  } else if (cmd=="to") {
    function(g) { igraph.es.to.igraphesseq(g, it) }
  } else if (cmd=="step") {
    igraph.es.next(NULL, it)
  } else {
    NULL
  }
}

"$.igraphesvector" <- function(it, cmd) {
  if (cmd=="end") {
    igraph.es.end(NULL, it)
  } else if (cmd=="e") {
    igraph.es.get.igraphesvector(NULL, it)
  } else if (cmd=="from") {
    function(g) { igraph.es.from.igraphesvector(g, it) }
  } else if (cmd=="to") {
    function(g) { igraph.es.to.igraphesvector(g, it) }
  } else if (cmd=="step") {
    igraph.es.next(NULL, it)
  } else {
    NULL
  }
}

#####
# internal helper functions

as.igraph.vs <- function(g, it) {
  if (class(it) %in% c("igraphvsseq", "igraphvsvector")) {
    it
  } else if (is.numeric(it)) {
    igraph.vs.vector(g, it)
  } else {
    stop("Cannot interpret as vertex set")
  }
}

as.igraph.es <- function(g, it) {
  if (class(it) %in% c("igraphesseq", "igraphesvector")) {
    it
  } else if (is.numeric(it)) {
    igraph.vs.vector(g, it)
  } else {
    stop("Cannot interpret as edge set")
  }
}
