
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
# Convert graphs to human readable forms
###################################################################

print.igraph <- function(x,
                         graph.attributes=igraph.par("print.graph.attributes"),
                         vertex.attributes=igraph.par("print.vertex.attributes"),
                         edge.attributes=igraph.par("print.edge.attributes"),
                         names=TRUE,
                         ...) {
  
  if (!is.igraph(x)) {
    stop("Not a graph object")
  }
  ec <- ecount(x)
  vc <- vcount(x)
  
  # From summary.graph
  cat("Vertices:", vc, "\n")
  cat("Edges:", ec, "\n")
  cat("Directed:", is.directed(x), "\n")

  # Graph attributes
  if (graph.attributes) {
    cat("Graph attributes:\n")
    list <- list.graph.attributes(x)
    sapply(list, function(n) {
      cat("  ", n, "=", get.graph.attribute(x, n), "\n") })
  }

  # Vertex attributes
  if (vertex.attributes) {
    cat("Vertex attributes:\n")
    list <- list.vertex.attributes(x)
    if (length(list) == 0) {
      cat("No vertex attributes\n")
    } else {
      if (vc==0 ||
          all(sapply(list, function(v) is.numeric(v) |
                     is.character(v) | is.logical(v)))) {
        ## create a table
        tab <- data.frame(v=paste(sep="", "[", seq(length=vc)-1, "]"), row.names="v")
        for (i in list) {
          tab[i] <- get.vertex.attribute(x, i)
        }
        print(tab)
      } else {
        for (i in 0:(vc-1)) {
          cat("  ", i, "  ")
          sapply(list, function(n) {
            cat(n, "=", get.vertex.attribute(x, n, i), "\t")})
          cat("\n")
        }
      }
    }
  }

  if (edge.attributes) {
    list <- list.edge.attributes(x)
  } else {
    list <- character()
  }
  
  arrow <- ifelse(is.directed(x), "->", "--")
  if (ec != 0) {
    if (!edge.attributes) {
      cat("Edges:\n")
    } else {
      cat("Edges and their attributes:\n")
    }
    el <- get.edgelist(x, names=names)
    if (ec==0 || 
        all(sapply(list, function(v) is.numeric(v) |
                   is.character(v) | is.logical(v)))) {
      ## create a table
      tab <- data.frame(e=paste(sep="", "[", seq(length=ec)-1, "]"), row.names="e")
      tab[" "] <- paste(el[,1], arrow, el[,2])
      for (i in list) {
        tab[i] <- get.edge.attribute(x, i)
      }
      print(tab)
    } else {
      i <- 0
      apply(el, 1, function(v) {
        cat(sep="", "[", i, "] ", v[1], " ", arrow, " ", v[2]);
        if (edge.attributes) {
          sapply(list, function(n) {
            cat("  ", n, "=", get.edge.attribute(x, n, i), "\t")})
        }
        cat("\n")
        i <<- i+1
      })
    }
  }
  
  invisible(x)
}

summary.igraph <- function(object, ...) {

  if (!is.igraph(object)) {
    stop("Not a graph object")
  }
  cat("Vertices:", vcount(object), "\n")
  cat("Edges:", ecount(object), "\n")
  cat("Directed:", is.directed(object), "\n")
  l <- list.graph.attributes(object)
  if (length(l)==0) {
    cat("No graph attributes.\n")
  } else {
    cat(sep="", "Graph attributes: ", paste(l, collapse=", "), ".\n")
  }
  l <- list.vertex.attributes(object)
  if (length(l)==0) {
    cat("No vertex attributes.\n")
  } else {
    cat(sep="", "Vertex attributes: ", paste(l, collapse=", "), ".\n")
  }
  l <- list.edge.attributes(object)
  if (length(l)==0) {
    cat("No edge attributes.\n")
  } else {
    cat(sep="", "Edge attributes: ", paste(l, collapse=", "), ".\n")
  }
  
  invisible(object)
}
