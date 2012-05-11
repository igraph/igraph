
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
#   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
#   02110-1301 USA
#
###################################################################

###################################################################
# Layouts
###################################################################

layout.random <- function(graph, params, dim=2) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (dim==2) {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    .Call("R_igraph_layout_random", graph,
          PACKAGE="igraph")
  } else if (dim==3) {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    .Call("R_igraph_layout_random_3d", graph,
          PACKAGE="igraph")
  } else {
    stop("Invalid `dim' value");
  }
}

layout.circle <- function(graph, params) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_layout_circle", graph,
        PACKAGE="igraph")
}

layout.sphere <- function(graph, params) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_layout_sphere", graph,
        PACKAGE="igraph")
}

layout.fruchterman.reingold <- function(graph, ..., dim=2,
                                        params=list()) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (length(params)==0) {
    params <- list(...)
  }

  if (dim==2) {
    fn <- "R_igraph_layout_fruchterman_reingold"
  } else if (dim==3 ){
    fn <- "R_igraph_layout_fruchterman_reingold_3d"
  } else {
    stop("Invalid `dim' argument");
  }
  
  vc <- vcount(graph)
  if (is.null(params$niter))     { params$niter      <- 500  }
  if (is.null(params$maxdelta))  { params$maxdelta   <- vc   }
  if (is.null(params$area))      { params$area       <- vc^2 }
  if (is.null(params$coolexp))   { params$coolexp    <- 1.5  }
  if (is.null(params$repulserad)){ params$repulserad <- params$area * vc }
  if (is.null(params$weights))   {
    params$weights <- NULL
  } else {
    params$weights <- as.numeric(params$weights)
  }
  if (!is.null(params$start)) {
    params$start <- structure(as.numeric(params$start), dim=dim(params$start))
  }
  if (!is.null(params$minx)) {
    params$minx <- as.double(params$minx)
  }
  if (!is.null(params$maxx)) {
    params$maxx <- as.double(params$maxx)
  }
  if (!is.null(params$miny)) {
    params$miny <- as.double(params$miny)
  }
  if (!is.null(params$maxy)) {
    params$maxy <- as.double(params$maxy)
  }
  if (!is.null(params$minz)) {
    params$minz <- as.double(params$minz)
  }
  if (!is.null(params$maxz)) {
    params$maxz <- as.double(params$maxz)
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call(fn, graph,
        as.double(params$niter), as.double(params$maxdelta),
        as.double(params$area), as.double(params$coolexp),
        as.double(params$repulserad), params$weights, params$start,
        params$minx, params$maxx, params$miny, params$maxy,
        params$minz, params$maxz,
        PACKAGE="igraph")
}

layout.fruchterman.reingold.grid <- function(graph, ...,
                                             params=list()) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (length(params)==0) {
    params <- list(...)
  }

  vc <- vcount(graph)
  if (is.null(params$niter))     { params$niter      <- 500  }
  if (is.null(params$maxdelta))  { params$maxdelta   <- vc   }
  if (is.null(params$area))      { params$area       <- vc^2 }
  if (is.null(params$coolexp))   { params$coolexp    <- 1.5  }
  if (is.null(params$repulserad)){ params$repulserad <- params$area * vc }
  if (is.null(params$cellsize))  { params$cellsize   <-
                                     (sqrt(sqrt(params$area))) }
  if (is.null(params$weights))   {
    params$weights <- NULL
  } else {
    params$weights <- as.numeric(params$weights)
  }
  if (!is.null(params$start)) {
    params$start <- structure(as.numeric(params$start), dim=dim(params$start))
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_layout_fruchterman_reingold_grid", graph,
        as.double(params$niter), as.double(params$maxdelta),
        as.double(params$area), as.double(params$coolexp),
        as.double(params$repulserad), as.double(params$cellsize),
        params$start, params$weights,
        PACKAGE="igraph")
}
  

# FROM SNA 0.5

layout.kamada.kawai<-function(graph, ..., dim=2,
                              params=list()) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (length(params)==0) {
    params <- list(...)
  }

  if (dim==2) {
    fn <- "R_igraph_layout_kamada_kawai"
  } else if (dim==3) {
    fn <- "R_igraph_layout_kamada_kawai_3d"
  } else {
    stop("Invalid `dim' parameter")
  }
  
  vc <- vcount(graph)
  if (is.null(params$niter))      { params$niter   <- 1000 }
  if (is.null(params$sigma))      { params$sigma   <- vc/4 }
  if (is.null(params$initemp))    { params$initemp <- 10   }
  if (is.null(params$coolexp))    { params$coolexp <- 0.99 }
  if (is.null(params$kkconst))    { params$kkconst <- vc^2 }
  if (is.null(params$fixz))       { params$fixz    <- FALSE}
  if (!is.null(params$start)) {
    params$start <- structure(as.numeric(params$start), dim=dim(params$start))
  }
  if (!is.null(params$minx)) {
    params$minx <- as.double(params$minx)
  }
  if (!is.null(params$maxx)) {
    params$maxx <- as.double(params$maxx)
  }
  if (!is.null(params$miny)) {
    params$miny <- as.double(params$miny)
  }
  if (!is.null(params$maxy)) {
    params$maxy <- as.double(params$maxy)
  }
  if (!is.null(params$minz)) {
    params$minz <- as.double(params$minz)
  }
  if (!is.null(params$maxz)) {
    params$maxz <- as.double(params$maxz)
  }
  if (params$fixz && dim==2) {
    warning("`fixz' works for 3D only, ignored.")
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call(fn, graph,
        as.double(params$niter), as.double(params$initemp),
        as.double(params$coolexp), as.double(params$kkconst),
        as.double(params$sigma), params$start, as.logical(params$fixz),
        params$minx, params$maxx, params$miny, params$maxy,
        params$minz, params$maxz,
        PACKAGE="igraph")
}

layout.graphopt <- function(graph, ..., 
                            params=list()) {
  
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (length(params)==0) {
    params <- list(...)
  }

  vc <- vcount(graph)
  if (is.null(params$niter))            { params$niter           <- 500   }
  if (is.null(params$charge))           { params$charge          <- 0.001 }
  if (is.null(params$mass))             { params$mass            <- 30    }
  if (is.null(params$spring.length))    { params$spring.length   <- 0     }
  if (is.null(params$spring.constant))  { params$spring.constant <- 1     }
  if (is.null(params$max.sa.movement))  { params$max.sa.movement <- 5     }
  if (!is.null(params$start)) {
    params$start <- structure(as.numeric(params$start), dim=dim(params$start))
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_layout_graphopt", graph,
        as.double(params$niter), as.double(params$charge),
        as.double(params$mass), as.double(params$spring.length),
        as.double(params$spring.constant), params$max.sa.movement,
        PACKAGE="igraph")
}

layout.lgl <- function(graph, ..., params=list()) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (length(params)==0) {
    params <- list(...)
  }

  vc <- vcount(graph)
  if (is.null(params$maxiter))   { params$maxiter    <- 150  }
  if (is.null(params$maxdelta))  { params$maxdelta   <- vc   }
  if (is.null(params$area))      { params$area       <- vc^2 }
  if (is.null(params$coolexp))   { params$coolexp    <- 1.5  }
  if (is.null(params$repulserad)){ params$repulserad <- params$area * vc }
  if (is.null(params$cellsize))  { params$cellsize   <-
                                     (sqrt(sqrt(params$area))) }
  if (is.null(params$root))      {
    params$root <- -1
  } else {
    params$root <- as.igraph.vs(graph, params$root)-1
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_layout_lgl", graph, as.double(params$maxiter),
        as.double(params$maxdelta), as.double(params$area),
        as.double(params$coolexp), as.double(params$repulserad),
        as.double(params$cellsize), params$root,
        PACKAGE="igraph")
}

layout.reingold.tilford <- function(graph, ..., params=list()) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (length(params)==0) {
    params <- list(...)
  }

  if (is.null(params$root))          { params$root       <- numeric()  }
  if (is.null(params$circular))      { params$circular   <- FALSE      }
  if (is.null(params$rootlevel))     { params$rootlevel  <- numeric()  }
  if (is.null(params$mode))          { params$mode       <- "out"      }
  if (is.null(params$flip.y))        { params$flip.y     <- TRUE       }
  params$mode <- tolower(params$mode)
  params$mode <- switch(params$mode, "out"=1, "in"=2, "all"=3, "total"=3)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_layout_reingold_tilford", graph,
               as.igraph.vs(graph, params$root)-1,
               as.double(params$mode), as.double(params$rootlevel),
               as.logical(params$circular),
               PACKAGE="igraph")
  if (params$flip.y) {
    res[,2] <- max(res[,2])-res[,2]
  }
  res
}

layout.merge <- function(graphs, layouts, method="dla") {

  if (!all(sapply(graphs, is.igraph))) {
    stop("Not a graph object")
  }
  if (method == "dla") {
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
    res <- .Call("R_igraph_layout_merge_dla",
                 graphs, layouts,
                 PACKAGE="igraph")
  } else {
    stop("Invalid `method'.")
  }
  res
}

# FROM SNA 0.5

symmetrize.mat <- function(mats,rule=c("weak", "strong", "lower", "upper")){
   rule <- igraph.match.arg(rule)
   
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

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
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

layout.norm <- function(layout, xmin=NULL, xmax=NULL, ymin=NULL, ymax=NULL,
                          zmin=NULL, zmax=NULL) {

  if (!is.matrix(layout)) {
    stop("`layout' not a matrix")
  }
  if (ncol(layout) != 2 && ncol(layout) != 3) {
    stop("`layout' should have 2 or three columns")
  }
  
  if (!is.null(xmin) && !is.null(xmax)) {
    layout[,1] <- .layout.norm.col(layout[,1], xmin, xmax)
  }

  if (!is.null(ymin) && !is.null(ymax)) {
    layout[,2] <- .layout.norm.col(layout[,2], ymin, ymax)
  }
  
  if (ncol(layout)==3 && !is.null(zmin) && !is.null(zmax)) {
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

layout.svd <- function(graph, d=shortest.paths(graph), ...)
  UseMethod("layout.svd", graph)

layout.svd.igraph <- function(graph, d=shortest.paths(graph), ...) {
    
    if (!is.igraph(graph)) {
      stop("Not a graph object")
    }

    clust <- clusters(graph)
    llist <- list()
    llen <- numeric()
    glist <- list()
    for(i in 1:length(clust$csize)){
        ind <- clust$membership==i
        
        if(length(which(ind))>=3){
            thisl <- svd(d[ind, ind], 2)[[2]]
            thisl[, 1] <- thisl[, 1]/dist(range(thisl[, 1]))
            thisl[, 2] <- thisl[, 2]/dist(range(thisl[, 2]))
            llist[i+1] <- list(thisl)
        }else if(length(which(ind))==2){
            llist[i+1] <- list(d[ind, ind])
        } else {
            llist[i+1] <- list(matrix(c(0, 0), nrow=1))
        }
        
        llen[i+1] <- length(which(ind))
        
        glist[i+1] <- list(induced.subgraph(graph, V(graph)[ind]))
    }
    
    ## merge them all:
    lmerged <- layout.merge(glist, llist)
    
    ## now reorder these rows to reflect original graph:
    l <- matrix(rep(NA, 2*vcount(graph)), ncol=2)
    l[order(clust$membership), ] <- lmerged
    return(l)
}

piecewise.layout <- function(graph, layout=layout.kamada.kawai, ...) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  
  V(graph)$id <- seq(vcount(graph))
  gl <- decompose.graph(graph)
  ll <- lapply(gl, layout, ...)
  
  l <- layout.merge(gl, ll)
  l[ unlist(sapply(gl, get.vertex.attribute, "id")), ] <- l[]
  l
}

layout.drl <- function(graph, use.seed = FALSE,
                       seed=matrix(runif(vcount(graph)*2), ncol=2),
                       options=igraph.drl.default,
                       weights=E(graph)$weight,
                       fixed=NULL,
                       dim=2)
{
    if (!is.igraph(graph)) {
        stop("Not a graph object")
    }
    if (dim != 2 && dim != 3) {
      stop("`dim' must be 2 or 3")
    }
    use.seed <- as.logical(use.seed)
    seed <- as.matrix(seed)
    options.tmp <- igraph.drl.default
    options.tmp[names(options)] <- options
    options <- options.tmp
    if (!is.null(weights)) {
      weights <- as.numeric(weights)
    }
    if (!is.null(fixed)) {
      fixed <- as.logical(fixed)
    }
    on.exit(.Call("R_igraph_finalizer", PACKAGE = "igraph"))
    if (dim==2) {
      res <- .Call("R_igraph_layout_drl", graph, seed, use.seed, options, 
                   weights, fixed, PACKAGE = "igraph")
    } else {
      res <- .Call("R_igraph_layout_drl_3d", graph, seed, use.seed, options, 
                   weights, fixed, PACKAGE = "igraph")
    }      
    res
}

igraph.drl.default <- list(edge.cut=32/40,
                           init.iterations=0,
                           init.temperature=2000,
                           init.attraction=10,
                           init.damping.mult=1.0,
                           liquid.iterations=200,
                           liquid.temperature=2000,
                           liquid.attraction=10,
                           liquid.damping.mult=1.0,
                           expansion.iterations=200,
                           expansion.temperature=2000,
                           expansion.attraction=2,
                           expansion.damping.mult=1.0,
                           cooldown.iterations=200,
                           cooldown.temperature=2000,
                           cooldown.attraction=1,
                           cooldown.damping.mult=.1,
                           crunch.iterations=50,
                           crunch.temperature=250,
                           crunch.attraction=1,
                           crunch.damping.mult=0.25,
                           simmer.iterations=100,
                           simmer.temperature=250,
                           simmer.attraction=.5,
                           simmer.damping.mult=0)

igraph.drl.coarsen <- list(edge.cut=32/40,
                           init.iterations=0,
                           init.temperature=2000,
                           init.attraction=10,
                           init.damping.mult=1.0,
                           liquid.iterations=200,
                           liquid.temperature=2000,
                           liquid.attraction=2,
                           liquid.damping.mult=1.0,
                           expansion.iterations=200,
                           expansion.temperature=2000,
                           expansion.attraction=10,
                           expansion.damping.mult=1.0,
                           cooldown.iterations=200,
                           cooldown.temperature=2000,
                           cooldown.attraction=1,
                           cooldown.damping.mult=.1,
                           crunch.iterations=50,
                           crunch.temperature=250,
                           crunch.attraction=1,
                           crunch.damping.mult=0.25,
                           simmer.iterations=100,
                           simmer.temperature=250,
                           simmer.attraction=.5,
                           simmer.damping.mult=0)

igraph.drl.coarsest <- list(edge.cut=32/40,
                            init.iterations=0,
                            init.temperature=2000,
                            init.attraction=10,
                            init.damping.mult=1.0,
                            liquid.iterations=200,
                            liquid.temperature=2000,
                            liquid.attraction=2,
                            liquid.damping.mult=1.0,
                            expansion.iterations=200,
                            expansion.temperature=2000,
                            expansion.attraction=10,
                            expansion.damping.mult=1.0,
                            cooldown.iterations=200,
                            cooldown.temperature=2000,
                            cooldown.attraction=1,
                            cooldown.damping.mult=.1,
                            crunch.iterations=200,
                            crunch.temperature=250,
                            crunch.attraction=1,
                            crunch.damping.mult=0.25,
                            simmer.iterations=100,
                            simmer.temperature=250,
                            simmer.attraction=.5,
                            simmer.damping.mult=0)

igraph.drl.refine <- list(edge.cut=32/40,
                          init.iterations=0,
                          init.temperature=50,
                          init.attraction=.5,
                          init.damping.mult=1.0,
                          liquid.iterations=0,
                          liquid.temperature=2000,
                          liquid.attraction=2,
                          liquid.damping.mult=1.0,
                          expansion.iterations=50,
                          expansion.temperature=500,
                          expansion.attraction=.1,
                          expansion.damping.mult=.25,
                          cooldown.iterations=50,
                          cooldown.temperature=250,
                          cooldown.attraction=1,
                          cooldown.damping.mult=.1,
                          crunch.iterations=50,
                          crunch.temperature=250,
                          crunch.attraction=1,
                          crunch.damping.mult=0.25,
                          simmer.iterations=0,
                          simmer.temperature=250,
                          simmer.attraction=.5,
                          simmer.damping.mult=0)

igraph.drl.final <- list(edge.cut=32/40,
                         init.iterations=0,
                         init.temperature=50,
                         init.attraction=.5,
                         init.damping.mult=0,
                         liquid.iterations=0,
                         liquid.temperature=2000,
                         liquid.attraction=2,
                         liquid.damping.mult=1.0,
                         expansion.iterations=50,
                         expansion.temperature=2000,
                         expansion.attraction=2,
                         expansion.damping.mult=1.0,
                         cooldown.iterations=50,
                         cooldown.temperature=200,
                         cooldown.attraction=1,
                         cooldown.damping.mult=.1,
                         crunch.iterations=50,
                         crunch.temperature=250,
                         crunch.attraction=1,
                         crunch.damping.mult=0.25,
                         simmer.iterations=25,
                         simmer.temperature=250,
                         simmer.attraction=.5,
                         simmer.damping.mult=0)

layout.auto <- function(graph, dim=2, ...) {

  ## 1. If there is a 'layout' graph attribute, we just use that.
  ## 2. Otherwise, if there are vertex attributes called 'x' and 'y',
  ##    we use those (and the 'z' vertex attribute as well, if present).
  ## 3. Otherwise, if the graph is connected and small (<100) we use
  ##    the Kamada-Kawai layout.
  ## 4. Otherwise if the graph is medium size (<1000) we use the
  ##    Fruchterman-Reingold layout.
  ## 5. Otherwise we use the DrL layout generator.
  
  if ("layout" %in% list.graph.attributes(graph)) {
    lay <- get.graph.attribute(graph, "layout")
    if (is.function(lay)) {
      lay(graph, ...)
    } else {
      lay
    }

  } else if ( all(c("x", "y") %in% list.vertex.attributes(graph)) ) {
    if ("z" %in% list.vertex.attributes(graph)) {
      cbind(V(graph)$x, V(graph)$y, V(graph)$z)
    } else {
      cbind(V(graph)$x, V(graph)$y)
    }

  } else if (is.connected(graph) && vcount(graph) < 100) {
    layout.kamada.kawai(graph, dim=dim, ...)

  } else if (vcount(graph) < 1000) {
    layout.fruchterman.reingold(graph, dim=dim, ...)

  } else {
    layout.drl(graph, dim=dim, ...)
  }
  
}

layout.sugiyama <- function(graph, layers=NULL, hgap=1, vgap=1,
                            maxiter=100, weights=NULL,
                            attributes=c("default", "all", "none")) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  if (!is.null(layers)) layers <- as.numeric(layers)-1
  hgap <- as.numeric(hgap)
  vgap <- as.numeric(vgap)
  maxiter <- as.integer(maxiter)
  if (is.null(weights) && "weight" %in% list.edge.attributes(graph)) { 
    weights <- E(graph)$weight 
  } 
  if (!is.null(weights) && any(!is.na(weights))) { 
    weights <- as.numeric(weights) 
  } else { 
    weights <- NULL 
  }
  attributes <- igraph.match.arg(attributes)
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_layout_sugiyama", graph, layers, hgap,
               vgap, maxiter, weights, PACKAGE="igraph")

  # Flip the y coordinates, more natural this way
  res$res[,2] <- max(res$res[,2]) - res$res[,2] + 1

  # Separate real and dummy vertices
  vc <- vcount(graph)
  res$layout <- res$res[seq_len(vc),]
  if (nrow(res$res)==vc) {
    res$layout.dummy <- matrix(nrow=0, ncol=2)
  } else {
    res$layout.dummy <- res$res[(vc+1):nrow(res$res),]
  }
  
  # Add some attributes to the extended graph
  E(res$extd_graph)$orig <- res$extd_to_orig_eids
  res$extd_to_orig_eids <- NULL

  res$extd_graph <- set.vertex.attribute(res$extd_graph, "dummy",
                                         value=c(rep(FALSE, vc),
                                           rep(TRUE, nrow(res$res)-vc)))

  res$extd_graph$layout <- rbind(res$layout, res$layout.dummy)

  if (attributes=="default" || attributes=="all") {
    if ("size" %in% list.vertex.attributes(graph)) {
      V(res$extd_graph)$size <- 0
      V(res$extd_graph)$size[ !V(res$extd_graph)$dummy ] <- V(graph)$size
    }
    if ("size2" %in% list.vertex.attributes(graph)) {
      V(res$extd_graph)$size2 <- 0
      V(res$extd_graph)$size2[ !V(res$extd_graph)$dummy ] <- V(graph)$size2
    }
    if ("shape" %in% list.vertex.attributes(graph)) {
      V(res$extd_graph)$shape <- "none"
      V(res$extd_graph)$shape[ !V(res$extd_graph)$dummy ] <- V(graph)$shape
    }
    if ("label" %in% list.vertex.attributes(graph)) {
      V(res$extd_graph)$label <- ""
      V(res$extd_graph)$label[ !V(res$extd_graph)$dummy ] <- V(graph)$label
    }
    if ("color" %in% list.vertex.attributes(graph)) {
      V(res$extd_graph)$color <- head(V(graph)$color, 1)
      V(res$extd_graph)$color[ !V(res$extd_graph)$dummy ] <- V(graph)$color
    }
    eetar <- get.edgelist(res$extd_graph, names=FALSE)[,2]
    E(res$extd_graph)$arrow.mode <- 0
    if ("arrow.mode" %in% list.edge.attributes(graph)) {
      E(res$extd_graph)$arrow.mode[ eetar <= vc ] <- E(graph)$arrow.mode
    } else {
      E(res$extd_graph)$arrow.mode[ eetar <= vc ] <- is.directed(graph) * 2
    }
    if ("arrow.size" %in% list.edge.attributes(graph)) {
      E(res$extd_graph)$arrow.size <- 0
      E(res$extd_graph)$arrow.size[ eetar <= vc ] <- E(graph)$arrow.size
    }
  }

  if (attributes=="all") {
    gatt <- setdiff(list.graph.attributes(graph), "layout")
    vatt <- setdiff(list.vertex.attributes(graph),
                    c("size", "size2", "shape", "label", "color"))
    eatt <- setdiff(list.edge.attributes(graph),
                    c("arrow.mode", "arrow.size"))
    for (ga in gatt) {
      res$extd_graph <- set.graph.attribute(res$extd_graph, ga,
                                            get.graph.attribute(graph, ga))
    }
    for (va in vatt) {
      notdummy <- which(!V(res$extd_graph)$dummy)
      res$extd_graph <- set.vertex.attribute(res$extd_graph, va,
                                             notdummy,
                                             get.vertex.attribute(graph, va))
    }
    for (ea in eatt) {
      eanew <- get.edge.attribute(graph, ea)[E(res$extd_graph)$orig]
      res$extd_graph <- set.edge.attribute(res$extd_graph, ea, value=eanew)
    }
  }
  
  res$res <- NULL
  res
}
