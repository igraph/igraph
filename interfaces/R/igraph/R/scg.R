
stochasticMatrix <- function(graph, norm = c("row", "col"), ...) {

  if (!is.igraph(graph)) {
    stop("Not a graph object")
  }

  norm <- igraph.match.arg(norm)
  
  graph <- get.adjacency(graph, ...)
  if (any(graph < 0)) { stop("'graph' must have non-negative entries") }

  if(norm == "row") {
    y <- rowSums(graph)
    if(any(y == 0)) { stop("zero out-degree in 'graph' not allowed") }
    return( graph/y )
  } else {
    y <- colSums(graph)
    if (any(y == 0)) { stop("zero in-degree in 'graph not allowed'") }
    return( t(t(graph) / y) )
  }
}

scgGrouping <- function(V, nt,
                         mtype=c("symmetric", "laplacian",
                           "stochastic"),
                         algo=c("optimum", "interv_km", "interv",
                           "exact_scg"),
                         p=NULL, maxiter=100) {

  V <- as.matrix(structure(as.double(V), dim=dim(V)))
  groups <- as.numeric(nt)

  mtype <- switch(igraph.match.arg(mtype), "symmetric"=1, 
                        "laplacian"=2, "stochastic"=3)
  algo <- switch(igraph.match.arg(algo), "optimum"=1,
                      "interv_km"=2, "interv"=3, "exact_scg"=4)
  if (!is.null(p)) p <- as.numeric(p)
  maxiter <- as.integer(maxiter)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_scg_grouping", V, as.integer(nt[1]),
               if (length(nt)==1) NULL else nt,
               mtype, algo, p, maxiter,
               PACKAGE="igraph")
  res
}

scgSemiProjectors <- function(groups,
                               mtype=c("symmetric", "laplacian",
                                 "stochastic"), p=NULL,
                               norm=c("row", "col"),
                               sparse=getIgraphOpt("sparsematrices")) {
  # Argument checks
  groups <- as.numeric(groups)-1
  mtype <- switch(igraph.match.arg(mtype), "symmetric"=1, 
  "laplacian"=2, "stochastic"=3)
  if (!is.null(p)) p <- as.numeric(p)
  norm <- switch(igraph.match.arg(norm), "row"=1, "col"=2)
  sparse <- as.logical(sparse)
  if (sparse && !require(Matrix)) { sparse <- FALSE }

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )
  # Function call
  res <- .Call("R_igraph_scg_semiprojectors", groups, mtype, p, norm,
               sparse,
               PACKAGE="igraph")

  if (sparse) {
    res$L <- igraph.i.spMatrix(res$L)
    res$R <- igraph.i.spMatrix(res$R)
  }
                
  res
}

scg <- function(X, ev, nt, groups=NULL, 
                mtype=c("symmetric", "laplacian", "stochastic"),
                algo=c("optimum", "interv_km", "interv",
                  "exact_scg"), norm=c("row", "col"),
                direction=c("default", "left", "right"),
                evec=NULL, p=NULL, use.arpack=FALSE, maxiter=300,
                sparse=getIgraphOpt("sparsematrices"),
                output=c("default", "matrix", "graph"), semproj=FALSE,
                epairs=FALSE, stat.prob=FALSE)
  UseMethod("scg")

scg.igraph <- function(X, ev, nt, groups=NULL,
                       mtype=c("symmetric", "laplacian", "stochastic"),
                       algo=c("optimum", "interv_km", "interv",
                         "exact_scg"), norm=c("row", "col"),
                       direction=c("default", "left", "right"),
                       evec=NULL, p=NULL, use.arpack=FALSE, maxiter=300,
                       sparse=getIgraphOpt("sparsematrices"),
                       output=c("default", "matrix", "graph"), semproj=FALSE,
                       epairs=FALSE, stat.prob=FALSE) {
  
  myscg(graph=X, matrix=NULL, sparsemat=NULL, ev=ev, nt=nt,
        groups=groups, mtype=mtype, algo=algo,
        norm=norm, direction=direction, evec=evec, p=p,
        use.arpack=use.arpack, maxiter=maxiter, sparse=sparse,
        output=output, semproj=semproj, epairs=epairs,
        stat.prob=stat.prob)
}

scg.matrix <- function(X, ev, nt, groups=NULL,
                       mtype=c("symmetric", "laplacian", "stochastic"),
                       algo=c("optimum", "interv_km", "interv",
                         "exact_scg"), norm=c("row", "col"),
                       direction=c("default", "left", "right"),
                       evec=NULL, p=NULL, use.arpack=FALSE, maxiter=300,
                       sparse=getIgraphOpt("sparsematrices"),
                       output=c("default", "matrix", "graph"), semproj=FALSE,
                       epairs=FALSE, stat.prob=FALSE) {
  
  myscg(graph=NULL, matrix=X, sparsemat=NULL, ev=ev, nt=nt,
        groups=groups, mtype=mtype, algo=algo,
        norm=norm, direction=direction, evec=evec, p=p, 
        use.arpack=use.arpack, maxiter=maxiter, sparse=sparse,
        output=output, semproj=semproj, epairs=epairs,
        stat.prob=stat.prob)
}

scg.Matrix <- function(X, ev, nt, groups=NULL,
                       mtype=c("symmetric", "laplacian", "stochastic"),
                       algo=c("optimum", "interv_km", "interv",
                         "exact_scg"), norm=c("row", "col"),
                       direction=c("default", "left", "right"),
                       evec=NULL, p=NULL, use.arpack=FALSE, maxiter=300,
                       sparse=getIgraphOpt("sparsematrices"),
                       output=c("default", "matrix", "graph"), semproj=FALSE,
                       epairs=FALSE, stat.prob=FALSE) {

  myscg(graph=NULL, matrix=NULL, sparsemat=X, ev=ev, nt=nt,
        groups=groups, mtype=mtype, algo=algo,
        norm=norm, direction=direction, evec=evec, p=p,
        use.arpack=use.arpack, maxiter=maxiter, sparse=sparse,
        output=output, semproj=semproj, epairs=epairs,
        stat.prob=stat.prob)
}

myscg <- function(graph, matrix, sparsemat, ev, nt, groups=NULL,
                  mtype=c("symmetric", "laplacian", "stochastic"),
                  algo=c("optimum", "interv_km", "interv",
                    "exact_scg"), norm=c("row", "col"),
                  direction=c("default", "left", "right"),
                  evec=NULL, p=NULL, use.arpack=FALSE, maxiter=300,
                  sparse=getIgraphOpt("sparsematrices"),
                  output=c("default", "matrix", "graph"), semproj=FALSE,
                  epairs=FALSE, stat.prob=FALSE) {

  ## Argument checks
  if (!is.null(graph))  { stopifnot(is.igraph(graph)) }
  if (!is.null(matrix)) { stopifnot(is.matrix(matrix)) }
  if (!is.null(sparsemat)) { stopifnot(inherits(sparsemat, "Matrix")) }

  if (!is.null(sparsemat)) { sparsemat <- as(sparsemat, "dgCMatrix") }
  ev <- as.numeric(as.integer(ev))
  nt <- as.numeric(as.integer(nt))
  if (!is.null(groups)) groups <- as.numeric(groups)
  mtype <- igraph.match.arg(mtype)
  algo <- switch(igraph.match.arg(algo), "optimum"=1,
                      "interv_km"=2, "interv"=3, "exact_scg"=4)
  if (!is.null(groups)) { storage.mode(groups) <- "double" }
  use.arpack <- as.logical(use.arpack)
  maxiter <- as.integer(maxiter)
  sparse <- as.logical(sparse)
  output <- switch(igraph.match.arg(output), "default"=1, "matrix"=2,
                   "graph"=3)
  semproj <- as.logical(semproj)
  epairs <- as.logical(epairs)

  on.exit( .Call("R_igraph_finalizer", PACKAGE="igraph") )

  if (mtype=="symmetric") {
    if (!is.null(evec)) { storage.mode(evec) <- "double" }
    res <- .Call("R_igraph_scg_adjacency", graph, matrix, sparsemat, ev,
                 nt, algo, evec, groups,
                 use.arpack, maxiter, sparse, output, semproj, epairs,
                 PACKAGE="igraph")

  } else if (mtype=="laplacian") {
    norm <- switch(igraph.match.arg(norm), "row"=1, "col"=2)
    if (!is.null(evec)) { storage.mode(evec) <- "complex" }
    direction <- switch(igraph.match.arg(direction), "default"=1, "left"=2,
                        "right"=3)
    res <- .Call("R_igraph_scg_laplacian", graph, matrix, sparsemat, ev,
                 nt, algo, norm, direction,
                 evec, groups, use.arpack, maxiter, sparse, output,
                 semproj, epairs,
                 PACKAGE="igraph")

  } else if (mtype=="stochastic") {
    norm <- switch(igraph.match.arg(norm), "row"=1, "col"=2)
    if (!is.null(evec)) { storage.mode(evec) <- "complex" }
    if (!is.null(p)) { storage.mode(p) <- "double" }
    stat.prob <- as.logical(stat.prob)
    res <- .Call("R_igraph_scg_stochastic", graph, matrix, sparsemat, ev,
                 nt, algo, norm, evec, groups, p, use.arpack,
                 maxiter, sparse, output, semproj, epairs, stat.prob,
                 PACKAGE="igraph")    
  }

  if (!is.null(res$Xt) &&
      class(res$Xt) == "igraph.tmp.sparse") {
    res$Xt <- igraph.i.spMatrix(res$Xt)
  }
  if (!is.null(res$L) && class(res$L) == "igraph.tmp.sparse") {
    res$L <- igraph.i.spMatrix(res$L)
  }
  if (!is.null(res$R) && class(res$R) == "igraph.tmp.sparse") {
    res$R <- igraph.i.spMatrix(res$R)
  }

  res
}
