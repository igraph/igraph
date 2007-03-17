
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

measure.dynamics.id <- function(graph, 
                                iterations=5,
                                time.window=NULL,
                                number=FALSE,
                                norm.method="old") {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  maxind <- max(degree(graph, mode="in"))

  st <- rep(1, vcount(graph))

  for (i in seq(along=numeric(iterations))) {

    # Standard deviation only at the last iteration
    if (i==iterations) {
      sd.real <- TRUE
    } else {
      sd.real <- FALSE
    }
    number.real <- number & sd.real
    
    if (i != 1) {
      if (norm.method=="old") {
        mes[[1]] <- mes[[1]]/mes[[1]][1]
      } else {
        mes[[1]] <- mes[[1]]/sum(mes[[1]])
      }
    }    

    if (is.null(time.window)) {
      mes <- .Call("R_igraph_measure_dynamics_id", graph,
                   as.numeric(st),
                   as.numeric(maxind), as.logical(sd.real),
                   as.logical(number.real),
                   PACKAGE="igraph")
    } else {
      mes <- .Call("R_igraph_measure_dynamics_idwindow", graph,
                   as.numeric(st),
                   as.numeric(time.window), as.numeric(maxind),
                   as.logical(sd.real),
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

  if (sd.real) {
    if (norm.method=="old") {
      mes[[2]] <- mes[[2]]/mes[[1]][1]
    } else {
      mes[[2]] <- mes[[2]]/sum(mes[[1]])
    }
  }

  if (norm.method=="old") {
    mes[[1]] <- mes[[1]]/mes[[1]][1]
  } else {
    mes[[1]] <- mes[[1]]/sum(mes[[1]])
  }

  if (number) {
    res <- list(akl=mes[[1]], st=st, sd=mes[[2]], no=mes[[3]])
  } else {
    res <- list(akl=mes[[1]], st=st, sd=mes[[2]])
  }
  
  res
}

measure.dynamics.idage <- function(graph, agebins=300,
                                   iterations=5,
                                   time.window=NULL,
                                   number=FALSE,
                                   norm.method="old") {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  maxind <- max(degree(graph, mode="in"))

  st <- rep(1, vcount(graph))
  sd <- TRUE

  if (!is.null(time.window)) {
    time.window <- as.numeric(time.window)
  }

  for (i in seq(along=numeric(iterations))) {

    # Standard deviation only at the last iteration
    if (sd && i==iterations) {
      sd.real <- TRUE
    } else {
      sd.real <- FALSE
    }
    number.real <- sd.real & number

    if (i != 1) {
      if (norm.method=="old") {
        mes[[1]] <- mes[[1]]/mes[[1]][1,1]
      } else {
        mes[[1]] <- mes[[1]]/sum(mes[[1]])
      }
    }    
    
    mes <- .Call("R_igraph_measure_dynamics_idage", graph,
                 as.numeric(st), as.numeric(agebins),
                 as.numeric(maxind), as.logical(sd.real),
                 as.logical(number.real),
                 time.window,
                 PACKAGE="igraph")

    mes[[1]][!is.finite(mes[[1]])] <- 0
    st <- .Call("R_igraph_measure_dynamics_idage_st", graph,
                mes[[1]], time.window,
                PACKAGE="igraph")
  }

#  print(mes[[1]][1,1])

  if (sd.real) {
    if (norm.method=="old") {
      mes[[2]] <- mes[[2]]/mes[[1]][1,1]
    } else {
      mes[[2]] <- mes[[2]]/sum(mes[[1]])
    }
  }
  if (norm.method=="old") {
    mes[[1]] <- mes[[1]]/mes[[1]][1,1]
  } else {
    mes[[1]] <- mes[[1]]/sum(mes[[1]])
  }

  if (number) {
    res <- list(akl=mes[[1]], st=st, sd=mes[[2]], no=mes[[3]])
  } else {
    res <- list(akl=mes[[1]], st=st, sd=mes[[2]])
  }
  
  res
}

measure.dynamics.citedcat.id.age <- function(graph, categories, agebins=300,
                                             iterations=5,
                                             norm=c(1,1,1)) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  maxind <- max(degree(graph, mode="in"))

  st <- rep(1, vcount(graph))
  sd <- TRUE                            # calculate sd?

  for (i in seq(along=numeric(iterations))) {

    # Standard deviation only at the last iteration
    if (sd && i==iterations) {
      sd.real <- TRUE
    } else {
      sd.real <- FALSE
    }

    if (i != 1) {
      mes[[1]] <- mes[[1]]/mes[[1]][norm[1],norm[2],norm[3]]
    }    

    mes <- .Call("R_igraph_measure_dynamics_citedcat_id_age", graph,
                 as.numeric(st), as.numeric(categories),
                 as.numeric(max(categories))+1,
                 as.numeric(agebins), as.numeric(maxind), as.logical(sd.real),
                 PACKAGE="igraph")
    
    mes[[1]][!is.finite(mes[[1]])] <- 0
    st <- .Call("R_igraph_measure_dynamics_citedcat_id_age_st", graph,
                mes[[1]], as.numeric(categories), 
                as.numeric(max(categories))+1,
                PACKAGE="igraph")
  }

#  print(mes[[1]][1,1])
  
  if (sd.real) {
    mes[[2]] <- mes[[2]]/mes[[1]][norm[1],norm[2],norm[3]]
  }
  mes[[1]] <- mes[[1]]/mes[[1]][norm[1],norm[2],norm[3]]

  res <- list(akl=mes[[1]], st=st, sd=mes[[2]])
  
  res
}  

measure.dynamics.citingcat.id.age <- function(graph, categories, agebins=300,
                                              iterations=5,
                                              norm=c(1,1,1)) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  maxind <- max(degree(graph, mode="in"))

  st <- rep(1, vcount(graph))
  sd <- TRUE

  for (i in seq(along=numeric(iterations))) {

    # Standard deviation only at the last iteration
    if (sd && i==iterations) {
      sd.real <- TRUE
    } else {
      sd.real <- FALSE
    }

    if (i != 1) {
      mes[[1]] <- mes[[1]]/mes[[1]][norm[1],norm[2],norm[3]]
    }    

    mes <- .Call("R_igraph_measure_dynamics_citingcat_id_age", graph,
                 as.numeric(st), as.numeric(categories),
                 as.numeric(max(categories))+1,
                 as.numeric(agebins), as.numeric(maxind), as.logical(sd.real),
                 PACKAGE="igraph")
    
    mes[[1]][!is.finite(mes[[1]])] <- 0
    st <- .Call("R_igraph_measure_dynamics_citingcat_id_age_st", graph,
                mes[[1]], as.numeric(categories), 
                as.numeric(max(categories))+1,
                PACKAGE="igraph")
  }

#  print(mes[[1]][1,1])
  
  if (sd.real) {
    mes[[2]] <- mes[[2]]/mes[[1]][norm[1],norm[2],norm[3]]
  }
  mes[[1]] <- mes[[1]]/mes[[1]][norm[1],norm[2],norm[3]]

  res <- list(akl=mes[[1]], st=st, sd=mes[[2]])  
  res
}  

## TODO: vtime and etime as attributes

measure.dynamics.d.d <- function(graph, vtime, etime,
                                 iterations=5) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  maxdeg <- max(degree(graph, mode="all"))
  events <- max(vtime, etime)+1
  sd <- TRUE
  
  st <- rep(1, events)

  for (i in seq(along=numeric(iterations))) {

#    print(i)
    
    # standard deviation only in the last iteration
    sd.real <- sd && i==iterations
    
    mes <- .Call("R_igraph_measure_dynamics_d_d", graph,
                 as.numeric(vtime), as.numeric(etime),
                 as.numeric(events), as.numeric(st),
                 as.numeric(maxdeg), as.logical(sd.real),
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

  res <- list(akl=mes[[1]], st=st, sd=mes[[2]])

  res
}
                 
measure.dynamics.lastcit <- function(graph, agebins, iterations=5,
                                     norm.method="old",
                                     number=FALSE) {

  if (!is.igraph(graph)) {
    stop("Not a graph object!")
  }

  st <- rep(1, vcount(graph))

  for (i in seq(along=numeric(iterations))) {

    ## Standard deviation only at the last iteration
    if (i==iterations) {
      sd.real <- TRUE
    } else {
      sd.real <- FALSE
    }
    number.real <- number & sd.real

    if (i != 1) {
      if (norm.method=="old") {
        mes[[1]] <- mes[[1]]/mes[[1]][1]
      } else {
        mes[[1]] <- mes[[1]]/sum(mes[[1]])
      }
    }

    mes <- .Call("R_igraph_measure_dynamics_lastcit", graph,
                 as.numeric(st), as.numeric(agebins), as.logical(sd.real),
                 as.logical(number.real),
                 PACKAGE="igraph")

    mes[[1]][!is.finite(mes[[1]])] <- 0

    st <- .Call("R_igraph_measure_dynamics_lastcit_st", graph, mes[[1]],
                PACKAGE="igraph")
  }

  if (sd.real) {
    if (norm.method=="old") {
      mes[[2]] <- mes[[2]]/mes[[1]][1]
    } else {
      mes[[2]] <- mes[[2]]/sum(mes[[1]])
    }
  }
  if (norm.method=="old") {
    mes[[1]] <- mes[[1]]/mes[[1]][1]
  } else {
    mes[[1]] <- mes[[1]]/sum(mes[[1]])
  }

  if (number) {
    res <- list(akl=mes[[1]], st=st, sd=mes[[2]], no=mes[[3]])
  } else {
    res <- list(akl=mes[[1]], st=st, sd=mes[[2]])
  }
  
  res
}

  
                                     
