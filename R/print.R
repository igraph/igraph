
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

print.graph <- function(graph) {

  # From summary.graph
  cat("Vertices:", vcount(graph), "\n")
  cat("Edges:", ecount(graph), "\n")
  cat("Directed:", is.directed(graph), "\n")
  cat("Type:", igraph.type(graph), "\n")

  arrow <- ifelse(is.directed(graph), "->", "--")
  if (ecount(graph) != 0) {
    cat("\nEdges:\n")
    for (i in 1:vcount(graph)) {
      neis <- neighbors(graph, i, "out")
      if (!is.directed(graph)) {
        no.loops <- sum(neis==i)
        neis <- c(neis[ neis > i ], rep(i, no.loops/2))
      }
      idx <- 1
      for (j in neis) {
        cat(sep="", "[", idx, "] ", i, " ", arrow, " ", j, "\n")
        idx <- idx + 1
      }
    }
  }
  
  invisible(graph)
}

summary.graph <- function(graph) {

  cat("Vertices:", vcount(graph), "\n")
  cat("Edges:", ecount(graph), "\n")
  cat("Directed:", is.directed(graph), "\n")
  cat("Type:", igraph.type(graph), "\n")
  
  invisible(graph)
}
