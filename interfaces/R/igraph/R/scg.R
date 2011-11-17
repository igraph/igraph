
stochasticMatrix <- function(graph, direction = c("row", "col"), ...) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  direction <- igraph.match.arg(direction)
  
  graph <- get.adjacency(graph, ...)
  if (any(graph < 0)) { stop("'graph' must have non-negative entries") }

  if(direction == "row") {
    y <- rowSums(graph)
    if(any(y == 0)) { stop("zero out-degree in 'graph' not allowed") }
    return( graph/y )
  } else {
    y <- colSums(graph)
    if (any(y == 0)) { stop("zero in-degree in 'graph not allowed'") }
    return( t(t(graph) / y) )
  }
}

scg.grouping <- function(V, groups,
                         matrix.type=c("symmetric", "laplacian",
                           "stochastic"),
                         algorithm=c("optimum", "interv_km", "interv",
                           "exact_scg"),
                         p=NULL, maxiter=100) {

  V <- structure(as.double(V), dim=dim(V))
  groups <- as.numeric(groups)

  matrix.type <- switch(igraph.match.arg(matrix.type), "symmetric"=1, 
                        "laplacian"=2, "stochastic"=3)
  algorithm <- switch(igraph.match.arg(algorithm), "optimum"=1,
                      "interv_km"=2, "interv"=3, "exact_scg"=4)
  if (!is.null(p)) p <- as.numeric(p)
  maxiter <- as.integer(maxiter)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_scg_grouping", V, as.integer(groups[1]),
               if (length(groups)==1) NULL else groups,
               matrix.type, algorithm, p, maxiter,
               PACKAGE="igraph")
  res
}

scg.semiProjectors <- function(groups,
                               matrix.type=c("symmetric", "laplacian",
                                 "stochastic"), p=NULL,
                               norm=c("row", "col"),
                               sparse=getIgraphOpt("sparsematrices")) {
  # Argument checks
  groups <- as.numeric(groups)-1
  matrix.type <- switch(igraph.match.arg(matrix.type), "symmetric"=1, 
  "laplacian"=2, "stochastic"=3)
  if (!is.null(p)) p <- as.numeric(p)
  norm <- switch(igraph.match.arg(norm), "row"=1, "col"=2)
  sparse <- as.logical(sparse)
  if (sparse && !require(Matrix)) { sparse <- FALSE }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_scg_semiprojectors", groups, matrix.type, p, norm,
               sparse,
               PACKAGE="igraph")

  if (sparse) {
    res$L <- igraph.i.spMatrix(res$L)
    res$R <- igraph.i.spMatrix(res$R)
  }
                
  res
}

scg <- function(X, ev, intervals, groups=NULL, 
                matrix.type=c("symmetric", "laplacian", "stochastic"),
                algorithm=c("optimum", "interv_km", "interv",
                  "exact_scg"), norm=c("row", "col"),
                direction=c("default", "left", "right"),
                evec=NULL, p=NULL, use.arpack=FALSE, maxiter=300)
  UseMethod("scg")

scg.igraph <- function(X, ev, intervals, groups=NULL,
                       matrix.type=c("symmetric", "laplacian", "stochastic"),
                       algorithm=c("optimum", "interv_km", "interv",
                         "exact_scg"), norm=c("row", "col"),
                       direction=c("default", "left", "right"),
                       evec=NULL, p=NULL, use.arpack=FALSE, maxiter=300) {
  
  myscg(graph=X, matrix=NULL, sparsemat=NULL, ev=ev, intervals=intervals,
        groups=groups, matrix.type=matrix.type, algorithm=algorithm,
        norm=norm, direction=direction, evec=evec, p=p,
        use.arpack=use.arpack, maxiter=maxiter)
}

scg.matrix <- function(X, ev, intervals, groups=NULL,
                       matrix.type=c("symmetric", "laplacian", "stochastic"),
                       algorithm=c("optimum", "interv_km", "interv",
                         "exact_scg"), norm=c("row", "col"),
                       direction=c("default", "left", "right"),
                       evec=NULL, p=NULL, use.arpack=FALSE, maxiter=300) {
  
  myscg(graph=NULL, matrix=X, sparsemat=NULL, ev=ev, intervals=intervals,
        groups=groups, matrix.type=matrix.type, algorithm=algorithm,
        norm=norm, direction=direction, evec=evec, p=p, 
        use.arpack=use.arpack, maxiter=maxiter)
}

scg.Matrix <- function(X, ev, intervals, groups=NULL,
                       matrix.type=c("symmetric", "laplacian", "stochastic"),
                       algorithm=c("optimum", "interv_km", "interv",
                         "exact_scg"), norm=c("row", "col"),
                       direction=c("default", "left", "right"),
                       evec=NULL, p=NULL, use.arpack=FALSE, maxiter=300) {

  myscg(graph=NULL, matrix=NULL, sparsemat=X, ev=ev, intervals=intervals,
        groups=groups, matrix.type=matrix.type, algorithm=algorithm,
        norm=norm, direction=direction, evec=evec, p=p, use.arpack=use.arpack,
        maxiter=maxiter)
}

myscg <- function(graph, matrix, sparsemat, ev, intervals, groups=NULL,
                  matrix.type=c("symmetric", "laplacian", "stochastic"),
                  algorithm=c("optimum", "interv_km", "interv",
                    "exact_scg"), norm=c("row", "col"),
                  direction=c("default", "left", "right"),
                  evec=NULL, p=NULL, use.arpack=FALSE, maxiter=300) {

  ## Argument checks
  if (!is.null(graph))  { stopifnot(is.igraph(graph)) }
  if (!is.null(matrix)) { stopifnot(is.matrix(matrix)) }
  if (!is.null(sparsemat)) { stopifnot(inherits(sparsemat, "Matrix")) }

  if (!is.null(sparsemat)) { sparsemat <- as(sparsemat, "dgCMatrix") }
  ev <- as.numeric(as.integer(ev))
  intervals <- as.numeric(as.integer(intervals))
  if (!is.null(groups)) groups <- as.numeric(groups)
  matrix.type <- igraph.match.arg(matrix.type)
  algorithm <- switch(igraph.match.arg(algorithm), "optimum"=1, "interv_km"=2,
                      "interv"=3, "exact_scg"=4)
  if (!is.null(groups)) { storage.mode(groups) <- "double" }
  use.arpack <- as.logical(use.arpack)
  maxiter <- as.integer(maxiter)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  if (matrix.type=="symmetric") {
    if (!is.null(evec)) { storage.mode(evec)   <- "double" }
    res <- .Call("R_igraph_scg_adjacency", graph, matrix, sparsemat, ev,
                 intervals=0L, intervals, algorithm, eval=NULL, evec, groups,
                 use.arpack, maxiter,
                 PACKAGE="igraph")
  } else if (matrix.type=="laplacian") {
    norm <- switch(igraph.match.arg(norm), "row"=1, "col"=2)
    if (!is.null(evec)) { storage.mode(evec)   <- "complex" }
    direction <- switch(igraph.match.arg(direction), "default"=1, "left"=2,
                        "right"=3)
    res <- .Call("R_igraph_scg_laplacian", graph, matrix, sparsemat, ev,
                 intervals=0L, intervals, algorithm, norm, direction,
                 eval=NULL, evec, groups, use.arpack, maxiter,
                 PACKAGE="igraph")
  } else if (matrix.type=="stochastic") {
    norm <- switch(igraph.match.arg(norm), "row"=1, "col"=2)
    if (!is.null(evec)) { storage.mode(evec)   <- "complex" }
    if (!is.null(p)) { storage.mode(p)      <- "double" }
    res <- .Call("R_igraph_scg_stochastic", graph, matrix, sparsemat, ev,
                 intervals=0L, intervals, algorithm, norm,
                 eval=NULL, evec, groups, p, use.arpack, maxiter,
                 PACKAGE="igraph")    
  }

  if (!is.null(res$scg_sparsemat)) {
    res$scg_sparsemat <- igraph.i.spMatrix(res$scg_sparsemat)
  }
  if (!is.null(res$Lsparse)) {
    res$Lsparse <- igraph.i.spMatrix(res$Lsparse)
  }
  if (!is.null(res$Rsparse)) {
    res$Rsparse <- igraph.i.spMatrix(res$Rsparse)
  }

  res
}
