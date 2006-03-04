
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

get.adjacency <- function(graph, type="both") {
  if (is.character(type)) {
    type <- switch(type, "upper"=0, "lower"=1, "both"=2)
  }
  
  .Call("R_igraph_get_adjacency", graph, as.numeric(type),
        PACKAGE="igraph")
}

get.edgelist <- function(graph) {
  matrix(.Call("R_igraph_get_edgelist", graph, TRUE,
               PACKAGE="igraph"), nc=2)
}
