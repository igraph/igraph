
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

get.adjacency <- function(graph) {

  vc <- vcount(graph)
  res <- matrix(0, vc, vc)

  if (vc != 0) {
    for (i in 1:vc) {
      neis <- neighbors(graph, i, "out")
      if (!is.directed(graph)) {
        loops <- sum(neis==i)
        neis <- neis[ neis!= i]
        neis <- c(neis, rep(i, loops/2))
      }
      for (n in neis) { res[i, n] <- res[i, n] + 1 }
    }
  }

  res
}

get.edgelist <- function(graph) {

  res <- matrix(0, ecount(graph), 2)
  vc <- vcount(graph)
  ix <- 1
  
  if (vc != 0) {
    for (i in 1:vc) {
      neis <- neighbors(graph, i, "out")
      if (!is.directed(graph)) {
        no.loops <- sum(neis==i)
        neis <- c(neis[ neis > i ], rep(i, no.loops/2))
      }
      res[seq(ix, length=length(neis)),1] <- i
      res[seq(ix, length=length(neis)),2] <- neis
      ix <- ix + length(neis)
    }
  }

  res
}
