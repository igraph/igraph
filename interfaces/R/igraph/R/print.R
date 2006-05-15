
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

print.igraph <- function(x, graph.attributes=FALSE,
                         vertex.attributes=FALSE, edge.attributes=FALSE,
                         ...) {
  
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
    if (vc != 0) {
      for (i in 0:(vc-1)) {
        cat("  ", i, "  ")
        sapply(list, function(n) {
          cat(n, "=", get.vertex.attribute(x, n, i), "\t")})
        cat("\n")
      }
    }
  }

  if (edge.attributes) {
    list <- list.edge.attributes(x)
  }
  
  arrow <- ifelse(is.directed(x), "->", "--")
  if (ec != 0) {
    if (!edge.attributes) {
      cat("\nEdges:\n")
    } else {
      cat("\nEdges and their attributes:\n")
    }
    i <- 0
    el <- get.edgelist(x)
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
  
  invisible(x)
}

summary.igraph <- function(object, ...) {

  cat("Vertices:", vcount(object), "\n")
  cat("Edges:", ecount(object), "\n")
  cat("Directed:", is.directed(object), "\n")
  
  invisible(object)
}
