
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

