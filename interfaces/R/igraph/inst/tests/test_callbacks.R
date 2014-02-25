
context("igraph callbacks in R")

test_that("igraph calls from callbacks are not allowed", {

  library(igraph)
  f <- function(graph, data, extra) { vcount(graph) }
  expect_that(graph.bfs(graph.ring(10), root=1, callback=f),
              throws_error("igraph callbacks cannot call igraph functions"))
})
