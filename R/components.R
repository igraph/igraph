
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
# Connected components, subgraphs, kinda
###################################################################

clusters <- function(graph, mode="weak") {

  if (is.directed(graph) && mode=="strong") {
    res <- .Call("REST_strong_components", igraph.c.interface, graph, 
                 PACKAGE="igraph")
  } else {
    res <- .Call("REST_clusters", igraph.c.interface, graph,
                 PACKAGE="igraph")
  }

  res
}

cluster.distribution <- function(graph, cumulative=FALSE, mul.size=FALSE,
                                 ...) {
  
  cs <- clusters(graph, ...)$csize;
  hi <- hist(cs, -1:max(cs), plot=FALSE)$intensities;
  if (mul.size) {
    hi <- hi*1:max(cs)
    hi <- hi/sum(hi)
  }
  if (!cumulative) {
    res <- hi
  } else {
    res <- rev(cumsum(rev(hi)));
  }
  
  res
}
