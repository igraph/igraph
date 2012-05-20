
#   IGraph R package
#   Copyright (C) 2006-2012  Gabor Csardi <csardi.gabor@gmail.com>
#   334 Harvard street, Cambridge, MA 02139 USA
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

graph.motifs <- function(graph, size=3, cut.prob=rep(0, size)) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  cut.prob <- as.numeric(cut.prob)
  if (length(cut.prob) != size) {
    cut.prob <- c(cut.prob[-length(cut.prob)],
                  rep(cut.prob[-length(cut.prob)], length(cut.prob)-1))
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_motifs_randesu", graph, as.integer(size),
        as.numeric(cut.prob),
        PACKAGE="igraph")
}

graph.motifs.no <- function(graph, size=3, cut.prob=rep(0, size)) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  cut.prob <- as.numeric(cut.prob)
  if (length(cut.prob) != size) {
    cut.prob <- c(cut.prob[-length(cut.prob)],
                  rep(cut.prob[-length(cut.prob)], length(cut.prob)-1))
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_motifs_randesu_no", graph, as.integer(size),
        as.numeric(cut.prob),
        PACKAGE="igraph")
}

graph.motifs.est <- function(graph, size=3, cut.prob=rep(0, size),
                             sample.size=vcount(graph)/10, sample=NULL) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  cut.prob <- as.numeric(cut.prob)
  if (length(cut.prob) != size) {
    cut.prob <- c(cut.prob[-length(cut.prob)],
                  rep(cut.prob[-length(cut.prob)], length(cut.prob)-1))
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_motifs_randesu_estimate", graph, as.integer(size),
        as.numeric(cut.prob), as.integer(sample.size), as.numeric(sample),
        PACKAGE="igraph")
}
  
