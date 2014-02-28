
#   IGraph R package
#   Copyright (C) 2003-2012  Gabor Csardi <csardi.gabor@gmail.com>
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

###################################################################
# Layouts
###################################################################

layout.random <- function(graph, dim=2) {
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

layout.circle <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_layout_circle", graph,
        PACKAGE="igraph")
}

layout.sphere <- function(graph) {
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_layout_sphere", graph,
        PACKAGE="igraph")
}  

layout.graphopt <- function(graph, start=NULL, niter=500, charge=0.001,
                            mass=30, sprint.length=0, spring.constant=1,
                            max.sa.movement=5) {
  
  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  if (!is.null(start)) {
    start <- structure(as.numeric(start), dim=dim(start))
  }
  niter <- as.double(niter)
  charge <- as.double(charge)
  mass <- as.double(mass)
  spring.length <- as.double(spring.length)
  spring.constant <- as.double(spring.constant)
  max.sa.movement <- as.double(max.sa.movement)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  .Call("R_igraph_layout_graphopt", graph, niter, charge, mass,
        spring.length, spring.constant, max.sa.movement, start,
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

layout.reingold.tilford <- function(graph, root=numeric(), circular=FALSE,
                                    rootlevel=numeric(), mode="out",
                                    flip.y=TRUE) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }
  root <- as.igraph.vs(graph, root)-1
  circular <- as.logical(circular)
  rootlevel <- as.double(rootlevel)
  mode <- switch(igraph.match.arg(mode), "out"=1, "in"=2, "all"=3,
                 "total"=3)
  flip.y <- as.logical(flip.y)
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  res <- .Call("R_igraph_layout_reingold_tilford", graph, root, mode,
               rootlevel, circular, PACKAGE="igraph")
  if (flip.y) { res[,2] <- max(res[,2])-res[,2] }
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
            llist[[i]] <- thisl
        }else if(length(which(ind))==2){
            llist[[i]] <- d[ind, ind]
        } else {
            llist[[i]] <- matrix(c(0, 0), nrow=1)
        }
        
        llen[i] <- length(which(ind))
        
        glist[[i]] <- induced.subgraph(graph, V(graph)[ind])
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
  ## 3. Otherwise, if the graph is small (<1000) we use
  ##    the Kamada-Kawai layout.
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

  } else if (vcount(graph) < 1000) {
    layout.kamada.kawai(graph, dim=dim, ...)

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

layout.mds <- function(graph, dist=NULL, dim=2,
                       options=igraph.arpack.default) {
  
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  if (!is.null(dist)) dist <- structure(as.double(dist), dim=dim(dist))
  dim <- as.integer(dim)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_layout_mds", graph, dist, dim,
        PACKAGE="igraph")

  res
}

layout.fruchterman.reingold <- function(graph, coords=NULL, dim=2,
                            niter=500, start.temp=sqrt(vcount(graph)),
                            grid=c("auto", "grid", "nogrid"), weights=NULL,
                            minx=NULL, maxx=NULL, miny=NULL, maxy=NULL,
                            minz=NULL, maxz=NULL,
                            coolexp, maxdelta, area, repulserad) {

                                        # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  if (!is.null(coords)) {
    coords <- as.matrix(structure(as.double(coords), dim=dim(coords)))
  }
  dim <- as.integer(dim)
  if (dim != 2L && dim != 3L) {
    stop("Dimension must be two or three")
  }
  niter <- as.integer(niter)
  start.temp <- as.numeric(start.temp)

  grid <- igraph.match.arg(grid)
  grid <- switch(grid, "grid"=0L, "nogrid"=1L, "auto"=2L)

  if (is.null(weights) && "weight" %in% list.edge.attributes(graph)) {
    weights <- E(graph)$weight
  }
  if (!is.null(weights) && any(!is.na(weights))) {
    weights <- as.numeric(weights)
  } else {
    weights <- NULL
  }
  if (!is.null(minx)) minx <- as.numeric(minx)
  if (!is.null(maxx)) maxx <- as.numeric(maxx)
  if (!is.null(miny)) miny <- as.numeric(miny)
  if (!is.null(maxy)) maxy <- as.numeric(maxy)
  if (!is.null(minz)) minz <- as.numeric(minz)
  if (!is.null(maxz)) maxz <- as.numeric(maxz)
  if (!missing(coolexp)) {
    warning("Argument `coolexp' is deprecated and has no effect")
  }
  if (!missing(maxdelta)) {
    warning("Argument `maxdelta' is deprecated and has no effect")
  }
  if (!missing(area)) {
    warning("Argument `area' is deprecated and has no effect")
  }
  if (!missing(repulserad)) {
    warning("Argument `repulserad' is deprecated and has no effect")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  if (dim==2) {
    res <- .Call("R_igraph_layout_fruchterman_reingold", graph, coords,
                 niter, start.temp, weights, minx, maxx, miny, maxy, grid,
                 PACKAGE="igraph")
  } else {
    res <- .Call("R_igraph_layout_fruchterman_reingold_3d", graph, coords,
                 niter, start.temp, weights, minx, maxx, miny, maxy,
                 minz, maxz, PACKAGE="igraph")
  }
  res
}

layout.kamada.kawai <- function(graph, coords=NULL, dim=2,
                                maxiter=50*vcount(graph),
                                epsilon=0.0, kkconst=vcount(graph),
                                weights=NULL, minx=NULL, maxx=NULL,
                                miny=NULL, maxy=NULL, minz=NULL, maxz=NULL,
                                niter, sigma, initemp, coolexp) {
  # Argument checks
  if (!is.igraph(graph)) { stop("Not a graph object") }
  if (!is.null(coords)) {
    coords <- as.matrix(structure(as.double(coords), dim=dim(coords)))
  }
  dim <- as.integer(dim)
  if (dim != 2L && dim != 3L) {
    stop("Dimension must be two or three")
  }

  maxiter <- as.integer(maxiter)
  epsilon <- as.numeric(epsilon)
  kkconst <- as.numeric(kkconst)
  if (is.null(weights) && "weight" %in% list.edge.attributes(graph)) { 
    weights <- E(graph)$weight 
  } 
  if (!is.null(weights) && any(!is.na(weights))) { 
    weights <- as.numeric(weights) 
  } else { 
    weights <- NULL 
  }
  if (!is.null(minx)) minx <- as.numeric(minx)
  if (!is.null(maxx)) maxx <- as.numeric(maxx)
  if (!is.null(miny)) miny <- as.numeric(miny)
  if (!is.null(maxy)) maxy <- as.numeric(maxy)
  if (!is.null(minz)) minz <- as.numeric(minz)
  if (!is.null(maxz)) maxz <- as.numeric(maxz)
  
  if (!missing(niter)) {
    warning("Argument `niter' is deprecated and has no effect")
  }
  if (!missing(sigma)) {
    warning("Argument `sigma' is deprecated and has no effect")
  }
  if (!missing(initemp)) {
    warning("Argument `initemp' is deprecated and has no effect")
  }
  if (!missing(coolexp)) {
    warning("Argument `coolexp' is deprecated and has no effect")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  if (dim == 2) {
    res <- .Call("R_igraph_layout_kamada_kawai", graph, coords, maxiter,
                 epsilon, kkconst, weights, minx, maxx, miny, maxy,
                 PACKAGE="igraph")
  } else {
    res <- .Call("R_igraph_layout_kamada_kawai_3d", graph, coords, maxiter,
                 epsilon, kkconst, weights, minx, maxx, miny, maxy, minz,
                 maxz, PACKAGE="igraph")
  }    

  res
}
