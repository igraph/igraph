
context("Bug 1019624")

test_that("weighted graph.adjacency works on integer matrices", {
  library(igraph)
  data <- matrix(c(0,0,0,2, 0,0,0,0, 0,0,0,2, 0,1,0,0), 4)
  g <- graph.adjacency(data, weighted=TRUE)
  expect_that(as.matrix(g[]), is_equivalent_to(data))
})


