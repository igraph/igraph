
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
# Connected components, subgraphs, kinda
###################################################################

no.clusters <- function(graph, mode=c("weak", "strong")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "weak"=1, "strong"=2)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_no_clusters", graph, as.numeric(mode),
        PACKAGE="igraph0")
}

cluster.distribution <- function(graph, cumulative=FALSE, mul.size=FALSE,
                                 ...) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  
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

is.connected <- function(graph, mode=c("weak", "strong")) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "weak"=1, "strong"=2)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_is_connected", graph, as.numeric(mode),
        PACKAGE="igraph0")
}

decompose.graph <- function(graph, mode=c("weak", "strong"), max.comps=NA,
                      min.vertices=0) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "weak"=1, "strong"=2)

  if (is.na(max.comps)) {
    max.comps=-1
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_decompose", graph, as.numeric(mode),
        as.numeric(max.comps), as.numeric(min.vertices),
        PACKAGE="igraph0"
        )
}
