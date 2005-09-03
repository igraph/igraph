
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


measure.dynamics.idage <- function(graph, agebins=300,
                                   iterations=5, sd=FALSE) {

  maxind <- max(degree(graph, mode="in"))

  st <- rep(1, vcount(graph))

  # Should we calculate the standard deviation?
  if (sd && iterations==0)
    { sd.real <- TRUE }
  else
    { sd.real <- FALSE }
  
  mes <- .Call("R_igraph_measure_dynamics_idage", graph,
               as.numeric(st), as.numeric(agebins),
               as.numeric(maxind), sd.real,
               PACKAGE="igraph")
  mes[[1]] <- mes[[1]]/mes[[1]][1,1]  
  
  for (i in seq(along=numeric(iterations))) {

    mes[[1]][!is.finite(mes[[1]])] <- 0
    st <- .Call("R_igraph_measure_dynamics_idage_st", graph,
                mes[[1]], 
                PACKAGE="igraph")

    # Should we also calculate the standard deviation?
    if (sd && i==tail(seq(along=numeric(iterations)),1) )
      { sd.real <- TRUE }
    else
      { sd.real <- FALSE }
    
    mes <- .Call("R_igraph_measure_dynamics_idage", graph,
                 as.numeric(st), as.numeric(agebins),
                 as.numeric(maxind), sd.real,
                 PACKAGE="igraph")
    if (sd.real) {
      mes[[2]] <- mes[[2]]/mes[[1]][1,1]
    }
#    mes[[3]][2:length(mes[[3]])] <- mes[[3]][2:length(mes[[3]])]/mes[[1]][1,1]
    mes[[1]] <- mes[[1]]/mes[[1]][1,1]
  }

  res <- list(akl=mes[[1]], st=st, sd=mes[[2]])
#  res$ize <- mes[[3]]
  res
}

