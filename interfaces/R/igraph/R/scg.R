
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
