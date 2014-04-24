
context("Bug 1073800")

test_that("Largest cliques is correct", {
  library(igraph)
  adj <- matrix(1, nrow=11, ncol=11) - diag(11)
  g <- graph.adjacency(adj)
  lc <- suppressWarnings(largest.cliques(g))
  expect_that(lc, equals(list(1:11)))
})
