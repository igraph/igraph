
#' Random walk on a graph
#'
#' Do a random walk. From the given start vertex, take the given number of
#' steps, choosing an edge from the actual vertex uniformly randomly. Edge
#' directions are observed in directed graphs (see the \code{mode} argument
#' as well). Multiple and loop edges are also observed.
#'
#' @param graph The input graph, might be undirected or directed.
#' @param start The start vertex.
#' @param steps The number of steps to make.
#' @param mode How to follow directed edges. \code{"out"} steps along the
#'   edge direction, \code{"in"} is opposite to that. \code{"all"} ignores
#'   edge directions. This argument is ignored for directed graphs.
#' @param stuck What to do if the random walk gets stuck. \code{"return"}
#'   returns the partial walk, \code{"error"} raises an error.
#' @return A vertex sequence containing the vertices along the walk.
#' @export
#' @examples
#' ## Stationary distribution of a Markov chain
#' g <- make_ring(10, directed = TRUE) %u%
#'   make_star(11, center = 11) + edge(11, 1)
#'
#' ec <- eigen_centrality(g, directed = TRUE)$vector
#' pg <- page_rank(g, damping = 0.999)$vector
#' w <- random_walk(g, start = 1, steps = 10000)
#'
#' ## These are similar, but not exactly the same
#' cor(table(w), ec)
#'
#' ## But these are (almost) the same
#' cor(table(w), pg)

random_walk <- function(graph, start, steps, mode = c("out", "in", "all"),
                        stuck = c("return", "error")) {
  ## Argument checks
  if (!is_igraph(graph)) stop("Not a graph object")
  start <- as.igraph.vs(graph, start)
  mode <- switch(igraph.match.arg(mode), "out" = 1, "in" = 2, "all" = 3,
                 "total" = 3)
  steps <- as.integer(steps)
  stuck <- switch(igraph.match.arg(stuck), "error" = 0L, "return" = 1L)

  on.exit( .Call("R_igraph_finalizer", PACKAGE = "igraph") )

  ## Function call
  res <- .Call("R_igraph_random_walk", graph, start - 1, mode, steps, stuck,
               PACKAGE="igraph")
  if (igraph_opt("return.vs.es")) {
    res <- create_vs(graph, res)
  }
  res
}
