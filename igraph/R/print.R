
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
    list <- g.a(x)
    sapply(list, function(n) { cat("  ", n, "=", g.a(x, n), "\n") })
  }

  # Vertex attributes
  if (vertex.attributes) {
    cat("Vertex attributes:\n")
    list <- v.a(x)
    if (vc != 0) {
      for (i in 0:(vc-1)) {
        cat("  ", i, "  ")
        sapply(list, function(n) { cat(n, "=", v.a(x, n, i), "\t")})
        cat("\n")
      }
    }
  }

  if (edge.attributes) {
    list <- e.a(x)
  }
  
  arrow <- ifelse(is.directed(x), "->", "--")
  if (ec != 0) {
    if (!edge.attributes) {
      cat("\nEdges:\n")
    } else {
      cat("\nEdges and their attributes:\n")
    }
    it <- ii.create(x, "eid")
    while (!ii.end(x, it)) {
      cat(sep="", "[", ii.get.edge(x, it), "] ",
          ii.get.from(x, it), " ", arrow, " ", ii.get.to(x, it))
      if (edge.attributes) {
        sapply(list, function(n)
               { cat("  ", n, "=", e.a(x, n, ii.get.edge(x,it)), "\t")})
      }
      cat("\n")
      it <- ii.next(x, it)
    }
  }
  
  invisible(x)
}

summary.igraph <- function(object, ...) {

  cat("Vertices:", vcount(object), "\n")
  cat("Edges:", ecount(object), "\n")
  cat("Directed:", is.directed(object), "\n")
  
  invisible(object)
}
