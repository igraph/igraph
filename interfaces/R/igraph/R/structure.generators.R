
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

graph <- function( edges, n=max(edges)+1, directed=TRUE ) {
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_create", as.numeric(edges), as.numeric(n),
        as.logical(directed),
        PACKAGE="igraph0")
}

graph.formula <- function(...) {
  mf <- as.list(match.call())[-1]

  ## Operators first
  f <- function(x) {
    if (is.call(x)) {
      return (list(as.character(x[[1]]), lapply(x[-1], f)))
    } else {
      return (NULL)
    }
  }
  ops <- unlist(lapply(mf, f))
  if (all(ops %in% c("-", ":"))) {
    directed <- FALSE
  } else if (all(ops %in% c("-", "+", ":"))) {
    directed <- TRUE
  } else {
    stop("Invalid operator in formula")
  }

  f <- function(x) {
    if (is.call(x)) {
      if (length(x)==3) {
        return( list(f(x[[2]]), op=as.character(x[[1]]), f(x[[3]])) )
      } else {
        return( list(op=as.character(x[[1]]), f(x[[2]])) )
      }
    } else {
      return( c(sym=as.character(x)) )
    }
  }
  
  ret <- lapply(mf, function(x) unlist(f(x)))

  v <- unique(unlist(lapply(ret, function(x) { x[ names(x)=="sym" ] })))

  ## Merge symbols for ":"
  ret <- lapply(ret, function(x) {
    res <- list()
    for (i in seq(along=x)) {
      if (x[i]==":" && names(x)[i]=="op") {
        ## SKIP
      } else if (i>1 && x[i-1]==":" && names(x)[i-1]=="op") {
        res[[length(res)]] <- c(res[[length(res)]], unname(x[i]))
      } else {
        res <- c(res, x[i])
      }
    }
    res
  })

  ## Ok, create the edges
  edges <- numeric()
  for (i in seq(along=ret)) {
    prev.sym <- character()
    lhead <- rhead <- character()
    for (j in seq(along=ret[[i]])) {
      act <- ret[[i]][[j]]
      if (names(ret[[i]])[j]=="op") {
        if (length(lhead)==0) {
          lhead <- rhead <- act
        } else {
          rhead <- act
        }
      } else if (names(ret[[i]])[j]=="sym") {
        for (ps in prev.sym) {
          for (ps2 in act) {
            if (lhead=="+") {
              edges <- c(edges, unname(c(ps2, ps)))
            }            
            if (!directed || rhead=="+") {
              edges <- c(edges, unname(c(ps, ps2)))
            }
          }
        }
        lhead <- rhead <- character()
        prev.sym <- act
      }
    }
  }

  ids <- seq(along=v)-1
  names(ids) <- v
  res <- graph( unname(ids[edges] ), n=length(v), directed=directed)
  res <- set.vertex.attribute(res, "name", value=v)
  res  
}

graph.adjacency.dense <- function(adjmatrix, mode=c("directed", "undirected", "max",
                                               "min", "upper", "lower", "plus"),
                                  weighted=NULL, diag=TRUE) {

  mode <- igraph.match.arg(mode)
  
  if (!diag) { diag(adjmatrix) <- 0 }
  
  if (!is.null(weighted)) {
    if (is.logical(weighted) && weighted) {
      weighted <- "weight"
    }
    if (!is.character(weighted)) {
      stop("invalid value supplied for `weighted' argument, please see docs.")
    }

    if (nrow(adjmatrix) != ncol(adjmatrix)) {
      stop("not a square matrix")
    }

    if (mode == "undirected") {
      if (!all(adjmatrix == t(adjmatrix))) {
        stop("Please supply a symmetric matrix if you want to create a weighted graph with mode=UNDIRECTED.")
      }
      adjmatrix[lower.tri(adjmatrix, diag=FALSE)] <- 0
    } else if (mode=="max") {
      adjmatrix <- pmax(adjmatrix, t(adjmatrix))
      adjmatrix[lower.tri(adjmatrix, diag=FALSE)] <- 0
    } else if (mode=="upper") {
      adjmatrix[lower.tri(adjmatrix, diag=FALSE)] <- 0
    } else if (mode=="lower") {
      adjmatrix[upper.tri(adjmatrix, diag=FALSE)] <- 0
    } else if (mode=="min") {
      adjmatrix <- pmin(adjmatrix, t(adjmatrix))
      adjmatrix[lower.tri(adjmatrix, diag=FALSE)] <- 0
    } else if (mode=="plus") {
      adjmatrix <- adjmatrix + t(adjmatrix)
      adjmatrix[lower.tri(adjmatrix, diag=FALSE)] <- 0
      diag(adjmatrix) <- diag(adjmatrix) / 2
    }
    
    no.edges <- sum(adjmatrix != 0)
    edges <- numeric(2*no.edges)
    weight <- numeric(no.edges)
    ptr <- 1
    if (no.edges == 0) {
      res <- graph.empty(directed=(mode==0))
      res <- set.edge.attribute(res, weighted, value=1)
      res
    } else {
      for (i in 1:nrow(adjmatrix)) {
        for (j in 1:ncol(adjmatrix)) {
          if (adjmatrix[i,j] != 0) {
            edges[2*ptr-1] <- i-1
            edges[2*ptr] <- j-1
            weight[ptr] <- adjmatrix[i,j]
            ptr <- ptr + 1
          }          
        }
      }
      res <- graph.empty(n=nrow(adjmatrix), directed=(mode=="directed"))
      weight <- list(weight)
      names(weight) <- weighted
      res <- add.edges(res, edges, attr=weight)
      res
    }
    
  } else {
    mode <- switch(mode,
                   "directed"=0,
                   "undirected"=1,
                   "max"=1,
                   "upper"=2,
                   "lower"=3,
                   "min"=4,
                   "plus"=5)
    
    adjmatrix <- as.matrix(adjmatrix)
    attrs <- attributes(adjmatrix)
    adjmatrix <- as.numeric(adjmatrix)
    attributes(adjmatrix) <- attrs
    
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
    res <- .Call("R_igraph_graph_adjacency", adjmatrix, as.numeric(mode),
                 PACKAGE="igraph0")
  }

  res
}

graph.adjacency.sparse <- function(adjmatrix, mode=c("directed", "undirected", "max",
                                                "min", "upper", "lower", "plus"),
                                   weighted=NULL, diag=TRUE) {

  require(Matrix)
  
  mode <- igraph.match.arg(mode)

  if (!is.null(weighted)) {
    if (is.logical(weighted) && weighted) {
      weighted <- "weight"
    }
    if (!is.character(weighted)) {
      stop("invalid value supplied for `weighted' argument, please see docs.")
    }
  }
  
  if (nrow(adjmatrix) != ncol(adjmatrix)) {
    stop("not a square matrix")
  }

  vc <- nrow(adjmatrix)

  ## to remove non-redundancies that can persist in a dgtMatrix
  if(is(adjmatrix, "dgTMatrix")) {
    adjmatrix = as(adjmatrix, "CsparseMatrix")
  }
  
  if (is.null(weighted) && mode=="undirected") { mode <- "max" }
  
  if (mode == "directed") {
    ## DIRECTED
    el <- selectMethod("summary", "sparseMatrix")(adjmatrix)
    if (!diag) { el <- el[ el[,1] != el[,2], ] }      
  } else if (mode == "undirected") {
    ## UNDIRECTED, must be symmetric if weighted
    if (!is.null(weighted) &&
        !selectMethod("isSymmetric", "sparseMatrix")(adjmatrix)) {
      stop("Please supply a symmetric matrix if you want to create a weighted graph with mode=UNDIRECTED.")
    }
    if (diag) {
      adjmatrix <- tril(adjmatrix)
    } else {
      adjmatrix <- tril(adjmatrix, -1)
    }      
    el <- selectMethod("summary", "sparseMatrix")(adjmatrix)
  } else if (mode=="max") {
    ## MAXIMUM
    el <- selectMethod("summary", "sparseMatrix")(adjmatrix)
    rm(adjmatrix)
    if (!diag) { el <- el[ el[,1] != el[,2], ] }
    el <- el[ el[,3] != 0, ]
    w <- el[,3]
    el <- el[,1:2]
    el <- cbind( pmin(el[,1],el[,2]), pmax(el[,1], el[,2]) )
    o <- order(el[,1], el[,2])
    el <- el[o,]
    w <- w[o]
    dd <- el[2:nrow(el),1] == el[1:(nrow(el)-1),1] &
          el[2:nrow(el),2] == el[1:(nrow(el)-1),2]
    dd <- which(dd)
    if (length(dd)>0) {
      mw <- pmax(w[dd], w[dd+1])
      w[dd] <- mw
      w[dd+1] <- mw
      el <- el[-dd,]
      w <- w[-dd]
    }
    el <- cbind(el, w)
    rm(w,dd)
  } else if (mode=="upper") {
    ## UPPER
    if (diag) {
      adjmatrix <- triu(adjmatrix)
    } else {
      adjmatrix <- triu(adjmatrix, 1)
    }
    el <- selectMethod("summary", "sparseMatrix")(adjmatrix)
    rm(adjmatrix)
    if (!diag) { el <- el[ el[,1] != el[,2], ] }      
  } else if (mode=="lower") {
    ## LOWER
    if (diag) {
      adjmatrix <- tril(adjmatrix)
    } else {
      adjmatrix <- tril(adjmatrix, -1)
    }
    el <- selectMethod("summary", "sparseMatrix")(adjmatrix)
    rm(adjmatrix)
    if (!diag) { el <- el[ el[,1] != el[,2], ] }      
  } else if (mode=="min") {
    ## MINIMUM
    adjmatrix <- sign(adjmatrix) * sign(Matrix::t(adjmatrix)) * adjmatrix
    el <- selectMethod("summary", "sparseMatrix")(adjmatrix)
    if (!diag) { el <- el[ el[,1] != el[,2], ] }
    el <- el[ el[,3] != 0, ]
    w <- el[,3]
    el <- el[,1:2]
    el <- cbind( pmin(el[,1],el[,2]), pmax(el[,1], el[,2]) )
    o <- order(el[,1], el[,2])
    el <- el[o,]
    w <- w[o]
    dd <- el[2:nrow(el),1] == el[1:(nrow(el)-1),1] &
          el[2:nrow(el),2] == el[1:(nrow(el)-1),2]
    dd <- which(dd)
    if (length(dd)>0) {
      mw <- pmin(w[dd], w[dd+1])
      w[dd] <- mw
      w[dd+1] <- mw
      el <- el[-dd,]
      w <- w[-dd]
    }
    el <- cbind(el, w)
    rm(w,dd)
  } else if (mode=="plus") {
    ## PLUS
    adjmatrix <- adjmatrix + Matrix::t(adjmatrix)
    if (diag) {
      adjmatrix <- tril(adjmatrix)
    } else {
      adjmatrix <- tril(adjmatrix, -1)
    }
    el <- selectMethod("summary", "sparseMatrix")(adjmatrix)
    if (diag) {
      loop <- el[,1] == el[,2]
      el[loop,3] <- el[loop,3] / 2
    }
    el <- el[ el[,3] != 0, ]
    rm(adjmatrix)
  }

  if (!is.null(weighted)) {
    res <- graph.empty(n=vc, directed=(mode=="directed"))
    weight <- list(el[,3])
    names(weight) <- weighted
    res <- add.edges(res, edges=t(as.matrix(el[,1:2]))-1, attr=weight)
  } else {
    edges <- unlist(apply(el, 1, function(x) rep(unname(x[1:2]), x[3])))
    res <- graph(n=vc, edges-1, directed=(mode=="directed"))
  }
  res
}

graph.adjacency <- function(adjmatrix, mode=c("directed", "undirected", "max",
                                         "min", "upper", "lower", "plus"),
                            weighted=NULL, diag=TRUE,
                            add.colnames=NULL, add.rownames=NA) {

  if (is(adjmatrix, "Matrix")) {
    res <- graph.adjacency.sparse(adjmatrix, mode=mode, weighted=weighted, diag=diag)
  } else {
    res <- graph.adjacency.dense(adjmatrix, mode=mode, weighted=weighted, diag=diag)
  }    
  
  ## Add columns and row names as attributes
  if (is.null(add.colnames)) {
    if (!is.null(colnames(adjmatrix))) {
      add.colnames <- "name"
    } else {
      add.colnames <- NA
    }
  } else if (!is.na(add.colnames)) {
    if (is.null(colnames(adjmatrix))) {
      warning("No column names to add")
      add.colnames <- NA
    }
  }
  
  if (is.null(add.rownames)) {
    if (!is.null(rownames(adjmatrix))) {
      add.rownames <- "name"
    } else {
      add.colnames <- NA
    }
  } else if (!is.na(add.rownames)) {
    if (is.null(rownames(adjmatrix))) {
      warning("No row names to add")
      add.rownames <- NA
    }
  }

  if (!is.na(add.rownames) && !is.na(add.colnames) &&
      add.rownames == add.colnames ) {
    warning("Same attribute for columns and rows, row names are ignored")
    add.rownames <- NA
  }

  if (!is.na(add.colnames)) {
    res <- set.vertex.attribute(res, add.colnames, value=colnames(adjmatrix))
  }
  if (!is.na(add.rownames)) {
    res <- set.vertex.attribute(res, add.rownames, value=rownames(adjmatrix))
  }

  res
}
  

graph.star <- function(n, mode=c("in", "out", "undirected"), center=0 ) {

  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=0, "in"=1, "undirected"=2)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_star", as.numeric(n), as.numeric(mode),
        as.numeric(center),
        PACKAGE="igraph0")
}

graph.full <- function(n, directed=FALSE, loops=FALSE) {
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_full", as.numeric(n), as.logical(directed),
        as.logical(loops),
        PACKAGE="igraph0")
}

###################################################################
# Lattices, every kind
###################################################################

graph.lattice <- function(dimvector=NULL,length=NULL, dim=NULL, nei=1,
                          directed=FALSE, mutual=FALSE, circular=FALSE, ...) {

##   # Check
##   if (is.null(dimvector) && (is.null(length) || is.null(dim))) {
##     stop("Either `length' and `dim' or 'dimvector' must be set. See docs.")
##   }
##   if (!is.null(length) && length < 1) {
##     stop("Invalid `length' argument, should be at least one")
##   }
##   if (!is.null(length) && dim < 1) {
##     stop("Invalid `dim' argument, should be at least one")
##   }
##   if (!is.null(length) && any(dimvector < 1)) {
##     stop("Invalid `dimvector', has negative or smaller than one elements")
##   }
##   if (mutual && !directed) {
##     warning("`mutual' specified for undirected graph, proceeding with multiplex edges...")
##   }
##   if (nei < 1) {
##     stop("`nei' should be at least one")
##   }
  
##   if (!is.null(length)) {
##     length <- as.numeric(length)
##     dim <- as.numeric(dim)
##     dimvector <- rep(length, times=dim)
##   } else {
##     dimvector <- as.numeric(dimvector)
##   }
##   nei <- as.numeric(nei)

##   n <- prod(dimvector)
##   res <- graph.empty(n=n, directed=directed, ...)
##   res <- add.edges(res, .Call("REST_create_lattice", dimvector, n,
##                               circular, mutual, PACKAGE="igraph0"))

##   # Connect also to local neighborhood
##   if (nei >= 2) {
##     neighbors <- lapply(1:length(res), function(a) get.neighborhood(res, a))
##     res <- add.edges(res, .Call("REST_connect_neighborhood", neighbors, nei,
##                                 mutual, PACKAGE="igraph0"))
##   }
  
##   res

  if (is.null(dimvector)) {
    dimvector <- rep(length, dim)
  }
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_lattice", as.numeric(dimvector), as.numeric(nei),
        as.logical(directed), as.logical(mutual),
        as.logical(circular),
        PACKAGE="igraph0")
}

graph.ring <- function(n, directed=FALSE, mutual=FALSE, circular=TRUE) {
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_ring", as.numeric(n), as.logical(directed),
        as.logical(mutual), as.logical(circular),
        PACKAGE="igraph0")
}

###################################################################
# Trees, regular
###################################################################

graph.tree <- function(n, children=2, mode=c("out", "in", "undirected")) {

  mode <- igraph.match.arg(mode)
  mode <- switch(mode, "out"=0, "in"=1, "undirected"=2);

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_tree", as.numeric(n), as.numeric(children),
        as.numeric(mode),
        PACKAGE="igraph0")
}

###################################################################
# The graph atlas
###################################################################

graph.atlas <- function(n) {

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_atlas", as.numeric(n),
        PACKAGE="igraph0")
}

###################################################################
# Create a graph from a data frame
###################################################################

graph.data.frame <- function(d, directed=TRUE, vertices=NULL) {

  if (ncol(d) < 2) {
    stop("the data frame should contain at least two columns")
  }

  ## Handle if some elements are 'NA'
  if (any(is.na(d[,1:2]))) {
    warning("In `d' `NA' elements were replaced with string \"NA\"")
    d[,1:2][ is.na(d[,1:2]) ] <- 'NA'
  }
  if (!is.null(vertices) && any(is.na(vertices[,1]))) {
    warning("In `vertices[,1]' `NA' elements were replaced with string \"NA\"")
    vertices[,1][is.na(vertices[,1])] <- 'NA'
  }    
  
  names <- unique( c(as.character(d[,1]), as.character(d[,2])) )
  if (!is.null(vertices)) {
    names2 <- names
    vertices <- as.data.frame(vertices)
    if (ncol(vertices) < 1) {
      stop("Vertex data frame contains no rows")
    }
    names <- as.character(vertices[,1])
    if (any(duplicated(names))) {
      stop("Duplicate vertex names")
    }
    if (any(! names2 %in% names)) {
      stop("Some vertex names in edge list are not listed in vertex data frame")
    }
  }
  ids <- seq(along=names)-1
  names(ids) <- names
    
  # create graph
  g <- graph.empty(n=0, directed=directed)

  # vertex attributes
  attrs <- list(name=names)
  if (!is.null(vertices)) {
    if (ncol(vertices) > 1) {
      for (i in 2:ncol(vertices)) {
        newval <- vertices[,i]
        if (class(newval) == "factor") {
          newval <- as.character(newval)
        }
        attrs[[ names(vertices)[i] ]] <- newval
      }
    }
  }

  # add vertices
  g <- add.vertices(g, length(ids), attr=attrs)
    
  # create edge list
  from <- as.character(d[,1])
  to <- as.character(d[,2])
  edges <- t(matrix(c(ids[from], ids[to]), ncol=2))
  
  # edge attributes
  attrs <- list()
  if (ncol(d) > 2) {
    for (i in 3:ncol(d)) {
      newval <- d[,i]
      if (class(newval) == "factor") {
        newval <- as.character(newval)
      }
      attrs[[ names(d)[i] ]] <- newval
    }
  }

  # add the edges
  g <- add.edges(g, edges, attr=attrs)
  g
}

graph.edgelist <- function(el, directed=TRUE) {

  if (!is.matrix(el) || ncol(el) != 2) {
    stop("graph.edgelist expects a matrix with two columns")
  }

  if (nrow(el) == 0) {
    res <- graph.empty(directed=directed)
  } else {  
    if (is.character(el)) {
      ## symbolic edge list
      names <- unique(as.character(t(el)))
      ids <- seq(names)-1
      names(ids) <- names
      res <- graph( unname(ids[t(el)]), directed=directed)
      rm(ids)
      V(res)$name <- names
    } else {
      ## normal edge list
      res <- graph( t(el), directed=directed )
    }
  }

  res
}

graph.extended.chordal.ring <- function(n, w) {
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_extended_chordal_ring", as.numeric(n),
        as.matrix(w),
        PACKAGE="igraph0")
}

line.graph <- function(graph) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_line_graph", graph,
        PACKAGE="igraph0")
}
  
  
graph.de.bruijn <- function(m, n) {

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_de_bruijn", as.numeric(m), as.numeric(n),
        PACKAGE="igraph0")
}

graph.kautz <- function(m, n) {

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_kautz", as.numeric(m), as.numeric(n),
        PACKAGE="igraph0")
}

graph.famous <- function(name) {

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  .Call("R_igraph_famous", as.character(name),
        PACKAGE="igraph0")
}

graph.full.bipartite <- function(n1, n2, directed=FALSE,
                                 mode=c("all", "out", "in")) {

  n1 <- as.numeric(n1)
  n2 <- as.numeric(n2)
  directed <- as.logical(directed)
  mode <- switch(igraph.match.arg(mode), "out"=1, "in"=2, "all"=3, "total"=3)  
  
  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  res <- .Call("R_igraph_full_bipartite", n1, n2, as.logical(directed), mode,
               PACKAGE="igraph0")
  set.vertex.attribute(res$graph, "type", value=res$types)
}

graph.bipartite <- function(types, edges, directed=FALSE) {

  types <- as.logical(types)
  edges <- as.numeric(edges)
  directed <- as.logical(directed)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
  res <- .Call("R_igraph_create_bipartite", types, edges, directed,
               PACKAGE="igraph0")
  set.vertex.attribute(res, "type", value=types)
}

graph.incidence.sparse <- function(incidence, directed, mode, multiple,
                                   weighted) {
  require(Matrix)

  n1 <- nrow(incidence)
  n2 <- ncol(incidence)
  el <- selectMethod("summary", signature=c(object="sparseMatrix"))(incidence)
  ## el <- summary(incidence)
  el[,2] <- el[,2] + n1
  el[,c(1,2)] <- el[,c(1,2)] - 1

  if (!is.null(weighted)) {

    if (is.logical(weighted) && weighted) {
      weighted <- "weight"
    }
    if (!is.character(weighted)) {
      stop("invalid value supplied for `weighted' argument, please see docs.")
    }

    if (!directed || mode==1) {
      ## nothing do to
    } else if (mode==2) {
      el[,1:2] <- el[,c(2,1)]
    } else if (mode==3) {
      el <- rbind(el, el[,c(2,1,3)])
    }

    res <- graph.empty(n=n1+n2, directed=directed)
    weight <- list(el[,3])
    names(weight) <- weighted
    res <- add.edges(res, edges=t(as.matrix(el[,1:2])), attr=weight)

  } else {

    if (multiple) {
      el[,3] <- ceiling(el[,3])
      el[,3][ el[,3] < 0 ] <- 0
    } else {
      el[,3] <- el[,3] != 0
    }

    if (!directed || mode==1) {
      ## nothing do to
    } else if (mode==2) {
      el[,1:2] <- el[,c(2,1)]
    } else if (mode==3) {
      el <- rbind(el, el[,c(2,1,3)])
    }
    
    edges <- unlist(apply(el, 1, function(x) rep(unname(x[1:2]), x[3])))
    res <- graph(n=n1+n2, edges, directed=directed)
  } 
    
  set.vertex.attribute(res, "type", value=c(rep(FALSE, n1), rep(TRUE, n2)))
}

graph.incidence.dense <- function(incidence, directed, mode, multiple,
                                  weighted) {
  
  if (!is.null(weighted)) {
    if (is.logical(weighted) && weighted) {
      weighted <- "weight"
    }
    if (!is.character(weighted)) {
      stop("invalid value supplied for `weighted' argument, please see docs.")
    }

    n1 <- nrow(incidence)
    n2 <- ncol(incidence)
    no.edges <- sum(incidence != 0)
    if (directed && mode==3) { no.edges <- no.edges * 2 }
    edges <- numeric(2*no.edges)
    weight <- numeric(no.edges)
    ptr <- 1
    for (i in seq_len(nrow(incidence))) {
      for (j in seq_len(ncol(incidence))) {
        if (incidence[i,j] != 0) {
          if (!directed || mode==1) {
            edges[2*ptr-1] <- i-1
            edges[2*ptr] <- n1+j-1
            weight[ptr] <- incidence[i,j]
            ptr <- ptr + 1
          } else if (mode==2) {
            edges[2*ptr-1] <- n1+j-1
            edges[2*ptr] <- i-1
            weight[ptr] <- incidence[i,j]
            ptr <- ptr + 1
          } else if (mode==3) {
            edges[2*ptr-1] <- i-1
            edges[2*ptr] <- n1+j-1
            weight[ptr] <- incidence[i,j]
            ptr <- ptr + 1
            edges[2*ptr-1] <- n1+j-1
            edges[2*ptr] <- i-1
          }
        }
      }
    }
    res <- graph.empty(n=n1+n2, directed=directed)
    weight <- list(weight)
    names(weight) <- weighted
    res <- add.edges(res, edges, attr=weight)
    res <- set.vertex.attribute(res, "type",
                                value=c(rep(FALSE, n1), rep(TRUE, n2)))
    
  } else {

    mode(incidence) <- "double"
    on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph0") )
    ## Function call
    res <- .Call("R_igraph_incidence", incidence, directed, mode, multiple,
                 PACKAGE="igraph0")
    res <- set.vertex.attribute(res$graph, "type", value=res$types)

  }

  res
}

graph.incidence <- function(incidence, directed=FALSE,
                            mode=c("all", "out", "in", "total"), 
                            multiple=FALSE, weighted=NULL,
                            add.names=NULL) {
  # Argument checks
  directed <- as.logical(directed)
  mode <- switch(igraph.match.arg(mode), "out"=1, "in"=2, "all"=3, "total"=3)
  multiple <- as.logical(multiple)

  if (is(incidence, "Matrix")) {
    res <- graph.incidence.sparse(incidence, directed=directed,
                                  mode=mode, multiple=multiple,
                                  weighted=weighted)
  } else {
    res <- graph.incidence.dense(incidence, directed=directed, mode=mode,
                                 multiple=multiple, weighted=weighted)
  }

  ## Add names
  if (is.null(add.names)) {
    if (!is.null(rownames(incidence)) && !is.null(colnames(incidence))) {
      add.names <- "name"
    } else {
      add.names <- NA
    }
  } else if (!is.na(add.names)) {
    if (is.null(rownames(incidence)) || is.null(colnames(incidence))) {
      warning("Cannot add row- and column names, at least one of them is missing")
      add.names <- NA
    }
  }
  if (!is.na(add.names)) {
    res <- set.vertex.attribute(res, add.names,
                                value=c(rownames(incidence), colnames(incidence)))
  }
  res
}

