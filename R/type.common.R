
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
# Common utilities
###################################################################

# ID generation

igraph.i.create.id <- function(id.length=20) {  
  chars <- c("a","b","c","d","e","f","g","h","i","j","k","l","m",
             "n","o","p","q","r","s","t","u","v","w","x","y","z",
             "A","B","C","D","E","F","G","H","I","J","K","L","M",
             "N","O","P","Q","R","S","T","U","V","W","X","Y","Z",
             "1","2","3","4","5","6","7","8","9","0")
  res <- paste(sample(chars, id.length, replace=TRUE), collapse="")

  res
} 

###################################################################
# Argument checks
###################################################################

add.edges.common <- function(graph, edges) {

  # Checks
  if (length(edges) %% 2 != 0) {
    stop("odd length of edges vector, should be even")
  }
  r <- range(edges)
  if (r[1]<1) {
    stop("Invalid (too small) vertex id's in edge set")
  }
  vc <- vcount(graph)
  if (r[2]>vc) {
    stop("Invalid (too big) vertex id's in edge set")
  }
  invisible(TRUE)
}

add.vertices.common <- function(graph, nv) {

  invisible(TRUE)
}

delete.edges.common <- function(graph, edges) {
  invisible(TRUE)
}
