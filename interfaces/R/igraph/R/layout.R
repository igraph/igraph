
#   IGraph R package
#   Copyright (C) 2003, 2004, 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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
# Layouts
###################################################################

layout.random <- function(graph, params) {
  .Call("R_igraph_layout_random", graph,
        PACKAGE="igraph")
}

layout.circle <- function(graph, params) {
  .Call("R_igraph_layout_circle", graph,
        PACKAGE="igraph")
}

layout.fruchterman.reingold <- function(graph, ..., params=list()) {

  if (length(params)==0) {
    params <- list(...)
  }

  vc <- vcount(graph)
  if (is.null(params$niter))     { params$niter      <- 500  }
  if (is.null(params$maxdelta))  { params$maxdelta   <- vc   }
  if (is.null(params$area))      { params$area       <- vc^2 }
  if (is.null(params$coolexp))   { params$coolexp    <- 1.5  }
  if (is.null(params$repulserad)){ params$repulserad <- params$area * vc }
  
  .Call("R_igraph_layout_fruchterman_reingold", graph,
        as.double(params$niter), as.double(params$maxdelta),
        as.double(params$area), as.double(params$coolexp),
        as.double(params$repulserad),
        PACKAGE="igraph")
}

# FROM SNA 0.5

layout.kamada.kawai<-function(graph, ..., params=list()) {

  if (length(params)==0) {
    params <- list(...)
  }

  vc <- vcount(graph)
  if (is.null(params$niter))      { params$niter   <- 1000 }
  if (is.null(params$sigma))      { params$sigma   <- vc/4 }
  if (is.null(params$initemp))    { params$initemp <- 10   }
  if (is.null(params$coolexp))    { params$coolexp <- 0.99 }
  if (is.null(params$kkconst))    { params$kkconst <- vc^2 }
  .Call("R_igraph_layout_kamada_kawai", graph,
        as.double(params$niter), as.double(params$initemp),
        as.double(params$coolexp), as.double(params$kkconst),
        as.double(params$sigma),
        PACKAGE="igraph")
}

layout.lgl <- function(graph, ..., params=list()) {

  if (length(params)==0) {
    params <- list(...)
  }

  vc <- vcount(graph)
  if (is.null(params$niter))     { params$maxiter    <- 500  }
  if (is.null(params$maxdelta))  { params$maxdelta   <- vc   }
  if (is.null(params$area))      { params$area       <- vc^2 }
  if (is.null(params$coolexp))   { params$coolexp    <- 1.5  }
  if (is.null(params$repulserad)){ params$repulserad <- params$area * vc }
  if (is.null(params$cellsize))  { params$cellsize   <- sqrt(params$area)/100 }
  
  .Call("R_igraph_layout_lgl", graph, as.double(params$maxiter),
        as.double(params$maxdelta), as.double(params$area),
        as.double(params$coolexp), as.double(params$repulserad),
        as.double(params$cellsize),
        PACKAGE="igraph")
}

# FROM SNA 0.5

symmetrize.mat <- function(mats,rule="weak"){
   #Build the input data structures
   if(length(dim(mats))>2){
      m<-dim(mats)[1]
      n<-dim(mats)[2]
      o<-dim(mats)[3]
      d<-mats
   }else{
      m<-1
      n<-dim(mats)[1]
      o<-dim(mats)[2]
      d<-array(dim=c(1,n,o))
      d[1,,]<-mats
   }
   #Apply the symmetry rule
   for(i in 1:m){
      if(rule=="upper"){
         temp<-d[i,,]
         for(j in 1:n)
            temp[j:n,j]<-temp[j,j:n]
         d[i,,]<-temp
      }else if(rule=="lower"){
         temp<-d[i,,]
         for(j in 1:n)
            temp[j,j:n]<-temp[j:n,j]
         d[i,,]<-temp
      }else if(rule=="weak"){
         d[i,,]<-matrix(as.numeric(d[i,,]|t(d[i,,])),nrow=n,ncol=o)
      }else if(rule=="strong"){
         d[i,,]<-matrix(as.numeric(d[i,,]&t(d[i,,])),nrow=n,ncol=o)
      }
   }
   #Return the symmetrized matrix
   if(m==1)
      out<-d[1,,]
   else
      out<-d
   out
}

# FROM SNA 0.5

layout.spring<-function(graph, ..., params=list()) {

  if (length(params)==0) {
    params <- list(...)
  }  

  if (is.null(params$mass))     { params$mass     <- 0.1 }
  if (is.null(params$equil))    { params$equil    <- 1 }
  if (is.null(params$k))        { params$k        <- 0.001 }
  if (is.null(params$repeqdis)) { params$repeqdis <- 0.1 }
  if (is.null(params$kfr))      { params$kfr      <- 0.01 }
  if (is.null(params$repulse))  { params$repulse  <- FALSE }

  #Create initial condidions
  vc <- vcount(graph)
  f.x <- rep(0,vc)       #Set initial x/y forces to zero
  f.y <- rep(0,vc)
  v.x <- rep(0,vc)       #Set initial x/y velocities to zero
  v.y <- rep(0,vc)
  tempa <- sample((0:(vc-1))/vc) #Set initial positions randomly on the circle
  x <- vc/(2*pi)*sin(2*pi*tempa)
  y <- vc/(2*pi)*cos(2*pi*tempa)
  ds <- symmetrize.mat(get.adjacency(graph))#Symmetrize/dichotomize the graph
  kfr <- params$kfr                     #Set initial friction level
  niter <- 1                            #Set the iteration counter
  #Simulate, with increasing friction, until motion stops    
  repeat{
    niter <- niter+1                    #Update the iteration counter
    dis <- as.matrix(dist(cbind(x,y)))  #Get inter-point distances
    #Get angles relative to the positive x direction
    theta <- acos(t(outer(x,x,"-"))/dis)*sign(t(outer(y,y,"-"))) 
    #Compute spring forces; note that we assume a base spring coefficient
    #of params$k units ("pseudo-Newtons/quasi-meter"?), with an equilibrium
    #extension of params$equil units for all springs
    f.x <- apply(ds*cos(theta)*params$k*(dis-params$equil),1,sum,na.rm=TRUE)
    f.y <- apply(ds*sin(theta)*params$k*(dis-params$equil),1,sum,na.rm=TRUE)
    #If node repulsion is active, add a force component for this
    #as well.  We employ an inverse cube law which is equal in power
    #to the attractive spring force at distance params$repeqdis
    if(params$repulse){
      f.x <- f.x-apply(cos(theta)*params$k/(dis/params$repeqdis)^3,1,
                       sum,na.rm=TRUE)
      f.y <- f.y-apply(sin(theta)*params$k/(dis/params$repeqdis)^3,1,
                       sum,na.rm=TRUE)
    }
    #Adjust the velocities (assume a mass of params$mass units); note that the
    #motion is roughly modeled on the sliding of flat objects across
    #a uniform surface (e.g., spring-connected cylinders across a table).
    #We assume that the coefficients of static and kinetic friction are
    #the same, which should only trouble you if you are under the 
    #delusion that this is a simulation rather than a graph drawing
    #exercise (in which case you should be upset that I'm not using
    #Runge-Kutta or the like!).
    v.x <- v.x+f.x/params$mass         #Add accumulated spring/repulsion forces
    v.y <- v.y+f.y/params$mass
    spd <- sqrt(v.x^2+v.y^2)     #Determine frictional forces
    fmag <- pmin(spd,kfr)  #We can't let friction _create_ motion!
    theta <- acos(v.x/spd)*sign(v.y)  #Calculate direction of motion
    f.x <- fmag*cos(theta)        #Decompose frictional forces
    f.y <- fmag*sin(theta)
    f.x[is.nan(f.x)] <- 0         #Correct for any 0/0 problems
    f.y[is.nan(f.y)] <- 0
    v.x <- v.x-f.x                #Apply frictional forces (opposing motion -
    v.y <- v.y-f.y                #note that mass falls out of equation)
    #Adjust the positions (yep, it's primitive linear updating time!)
    x <- x+v.x
    y <- y+v.y
    #Check for cessation of motion, and increase friction
    mdist <- mean(dis)
    if(all(v.x<mdist*1e-5)&&all(v.y<mdist*1e-5))
      break
    else
      kfr <- params$kfr*exp(0.1*niter)
  }
  #Return the result
  cbind(x,y)
}

i.layout.norm <- function(layout, xmin=NULL, xmax=NULL, ymin=NULL, ymax=NULL,
                          zmin=NULL, zmax=NUL) {

  if (!is.null(xmin) && !is.null(xmax)) {
    layout[,1] <- .layout.norm.col(layout[,1], xmin, xmax)
  }

  if (!is.null(ymin) && !is.null(ymax)) {
    layout[,2] <- .layout.norm.col(layout[,2], ymin, ymax)
  }
  
  if (!is.null(zmin) && !is.null(zmax)) {
    layout[,3] <- .layout.norm.col(layout[,3], zmin, zmax)
  }

  layout
}

.layout.norm.col <- function(v, min, max) {

  vr <- range(v)
  if (vr[1]==vr[2]) {
    fac <- 1
  } else {
    fac <- (max-min)/(vr[2]-vr[1])
  }

  (v-vr[1]) * fac + min
}
