
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

measure.dynamics.id <- function(graph, start.vertex=0,
                                iterations=5, significance=0,
                                estind=NULL, estage=NULL,
                                number=FALSE, time.window=NULL) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  maxind <- max(degree(graph, mode="in"))

  st <- rep(1, vcount(graph))
  sd <- (significance != 0)

  if (sd) {
    warning("Error calculation not implemented yet")
  }

  for (i in seq(along=numeric(iterations))) {

    # Standard deviation only at the last iteration
    if (sd && i==iterations) {
      sd.real <- significance
    } else {
      sd.real <- 0.0
    }

    # Number of estimates also
    if (number && i==iterations) {
      number.real <- number
    } else {
      number.real <- FALSE
    }

    if (i != 1) {
      mes[[1]] <- mes[[1]]/mes[[1]][1]
    }    

    if (is.null(time.window)) {
      mes <- .Call("R_igraph_measure_dynamics_id", graph,
                   as.numeric(start.vertex),
                   as.numeric(st),
                   as.numeric(maxind), as.numeric(sd.real),
                   as.logical(number.real),
                   PACKAGE="igraph")
    } else {
      mes <- .Call("R_igraph_measure_dynamics_idwindow", graph,
                   as.numeric(start.vertex), as.numeric(st),
                   as.numeric(time.window), as.numeric(maxind),
                   as.numeric(sd.real), as.logical(number.real),
                   PACKAGE="igraph")
    }

    mes[[1]][!is.finite(mes[[1]])] <- 0
    if (is.null(time.window)) {
      st <- .Call("R_igraph_measure_dynamics_id_st", graph,
                  mes[[1]], 
                  PACKAGE="igraph")
    } else {
      st <- .Call("R_igraph_measure_dynamics_idwindow_st", graph,
                  mes[[1]], as.numeric(time.window),
                  PACKAGE="igraph")
    }
  }

#  print(mes[[1]][1])
  
  if (sd.real != 0) {
    mes[[2]] <- mes[[2]]/mes[[1]][1]
    mes[[3]] <- mes[[3]]/mes[[1]][1]
  }
  if (!is.null(estind) && !is.null(estage)) {
    mes[[4]] <- mes[[4]]/mes[[1]][1]
##    mes[[4]] <- mes[[4]][-length(mes[[4]])]
  }
  mes[[1]] <- mes[[1]]/mes[[1]][1]

  res <- list(akl=mes[[1]], st=st, sd=mes[[2]],
              error=mes[[3]], no=mes[[4]], est=NULL)
  
  res
}

measure.dynamics.idage <- function(graph, start.vertex=0, agebins=300,
                                   iterations=5, significance=0,
                                   estind=NULL, estage=NULL,
                                   number=FALSE, time.window=NULL) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  maxind <- max(degree(graph, mode="in"))

  st <- rep(1, vcount(graph))
  sd <- (significance != 0)

  if (!is.null(time.window)) {
    time.window <- as.numeric(time.window)
  }

  for (i in seq(along=numeric(iterations))) {

    # Standard deviation only at the last iteration
    if (sd && i==iterations) {
      sd.real <- significance
    } else {
      sd.real <- 0.0
    }

    # Number of estimates also
    if (number && i==iterations) {
      number.real <- number
    } else {
      number.real <- FALSE
    }

    if (i != 1) {
      mes[[1]] <- mes[[1]]/mes[[1]][1,1]
    }    
    
    if (i==iterations && !is.null(estind) && !is.null(estage)) {
      mes <- .Call("R_igraph_measure_dynamics_idage_debug", graph,
                   as.numeric(st), as.numeric(agebins),
                   as.numeric(maxind), as.numeric(sd.real), as.numeric(estind),
                   as.numeric(estage), as.logical(number.real),
                   PACKAGE="igraph")
    } else {
      mes <- .Call("R_igraph_measure_dynamics_idage", graph,
                   as.numeric(start.vertex),
                   as.numeric(st), as.numeric(agebins),
                   as.numeric(maxind), as.numeric(sd.real),
                   as.logical(number.real), time.window,
                   PACKAGE="igraph")
    }

    mes[[1]][!is.finite(mes[[1]])] <- 0
    st <- .Call("R_igraph_measure_dynamics_idage_st", graph,
                mes[[1]], time.window,
                PACKAGE="igraph")
  }

#  print(mes[[1]][1,1])
  
  if (sd.real != 0) {
    mes[[2]] <- mes[[2]]/mes[[1]][1,1]
    mes[[3]] <- mes[[3]]/mes[[1]][1,1]
  }
  if (!is.null(estind) && !is.null(estage)) {
    mes[[4]] <- mes[[4]]/mes[[1]][1,1]
##    mes[[4]] <- mes[[4]][-length(mes[[4]])]
  }
  mes[[1]] <- mes[[1]]/mes[[1]][1,1]

  if (!is.null(estind) && !is.null(estage)) {
    res <- list(akl=mes[[1]], st=st, sd=mes[[2]],
                error=mes[[3]], no=mes[[4]], est=mes[[5]])
  } else {
    res <- list(akl=mes[[1]], st=st, sd=mes[[2]],
                error=mes[[3]], no=mes[[4]], est=NULL)
  }
  
  res
}

measure.dynamics.citedcat.id.age <- function(graph, categories, agebins=300,
                                             iterations=5, significance=0,
                                             number=FALSE, norm=c(1,1,1)) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  maxind <- max(degree(graph, mode="in"))

  st <- rep(1, vcount(graph))
  sd <- (significance != 0)

  for (i in seq(along=numeric(iterations))) {

    # Standard deviation only at the last iteration
    if (sd && i==iterations) {
      sd.real <- significance
    } else {
      sd.real <- 0.0
    }

    # Number of estimates also
    if (number && i==iterations) {
      number.real <- number
    } else {
      number.real <- FALSE
    }

    if (i != 1) {
      mes[[1]] <- mes[[1]]/mes[[1]][norm[1],norm[2],norm[3]]
    }    

    mes <- .Call("R_igraph_measure_dynamics_citedcat_id_age", graph,
                 as.numeric(0),
                 as.numeric(st), as.numeric(categories),
                 as.numeric(max(categories))+1,
                 as.numeric(agebins), as.numeric(maxind), as.numeric(sd.real),
                 as.logical(number.real),
                 PACKAGE="igraph")
    
    mes[[1]][!is.finite(mes[[1]])] <- 0
    st <- .Call("R_igraph_measure_dynamics_citedcat_id_age_st", graph,
                mes[[1]], as.numeric(categories), 
                as.numeric(max(categories))+1,
                PACKAGE="igraph")
  }

#  print(mes[[1]][1,1])
  
  if (sd.real != 0) {
    mes[[2]] <- mes[[2]]/mes[[1]][norm[1],norm[2],norm[3]]
    mes[[3]] <- mes[[3]]/mes[[1]][norm[1],norm[2],norm[3]]
  }
  mes[[1]] <- mes[[1]]/mes[[1]][norm[1],norm[2],norm[3]]

  res <- list(akl=mes[[1]], st=st, sd=mes[[2]],
              error=mes[[3]], no=mes[[4]], est=NULL)
  
  res
}  

measure.dynamics.citingcat.id.age <- function(graph, categories, agebins=300,
                                              iterations=5, significance=0,
                                              number=FALSE, norm=c(1,1,1)) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  maxind <- max(degree(graph, mode="in"))

  st <- rep(1, vcount(graph))
  sd <- (significance != 0)

  for (i in seq(along=numeric(iterations))) {

    # Standard deviation only at the last iteration
    if (sd && i==iterations) {
      sd.real <- significance
    } else {
      sd.real <- 0.0
    }

    # Number of estimates also
    if (number && i==iterations) {
      number.real <- number
    } else {
      number.real <- FALSE
    }

    if (i != 1) {
      mes[[1]] <- mes[[1]]/mes[[1]][norm[1],norm[2],norm[3]]
    }    

    mes <- .Call("R_igraph_measure_dynamics_citingcat_id_age", graph,
                 as.numeric(0),
                 as.numeric(st), as.numeric(categories),
                 as.numeric(max(categories))+1,
                 as.numeric(agebins), as.numeric(maxind), as.numeric(sd.real),
                 as.logical(number.real),
                 PACKAGE="igraph")
    
    mes[[1]][!is.finite(mes[[1]])] <- 0
    st <- .Call("R_igraph_measure_dynamics_citingcat_id_age_st", graph,
                mes[[1]], as.numeric(categories), 
                as.numeric(max(categories))+1,
                PACKAGE="igraph")
  }

#  print(mes[[1]][1,1])
  
  if (sd.real != 0) {
    mes[[2]] <- mes[[2]]/mes[[1]][norm[1],norm[2],norm[3]]
    mes[[3]] <- mes[[3]]/mes[[1]][norm[1],norm[2],norm[3]]
  }
  mes[[1]] <- mes[[1]]/mes[[1]][norm[1],norm[2],norm[3]]

  res <- list(akl=mes[[1]], st=st, sd=mes[[2]],
              error=mes[[3]], no=mes[[4]], est=NULL)
  
  res
}  

## TODO: vtime and etime as attributes

measure.dynamics.d.d <- function(graph, vtime, etime,
                                 iterations=5, sd=TRUE, no=FALSE) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  maxdeg <- max(degree(graph, mode="all"))
  events <- max(vtime, etime)+1

  st <- rep(1, events)

  for (i in seq(along=numeric(iterations))) {

#    print(i)
    
    # standard deviation only in the last iteration
    sd.real <- sd && i==iterations
    no.real <- no && i==iterations
    
    mes <- .Call("R_igraph_measure_dynamics_d_d", graph,
                 as.numeric(vtime), as.numeric(etime),
                 as.numeric(events), as.numeric(st),
                 as.numeric(maxdeg), as.logical(sd.real), as.logical(no.real),
                 PACKAGE="igraph")

    mes[[1]][!is.finite(mes[[1]])] <- 0 ## Hmmmm
    
    mes[[1]] <- mes[[1]]/mes[[1]][1]
    if (sd.real) {
      mes[[2]] <- mes[[2]]/mes[[1]][1]
    }

    st <- .Call("R_igraph_measure_dynamics_d_d_st", graph,
                as.numeric(vtime), as.numeric(etime),
                mes[[1]], as.numeric(events), as.numeric(maxdeg),
                PACKAGE="igraph")
  }

  res <- list(akl=mes[[1]], st=st, sd=mes[[2]], no=mes[[3]])

  res
}
                 
                 
