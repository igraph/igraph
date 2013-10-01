
context("graph.atlas")

test_that("graph.atlas works", {
  library(igraph)
  g124 <- graph.atlas(124)
  expect_that(graph.isomorphic(g124, graph(c(1,2,2,3,3,4,4,5,1,5,1,3,2,6),
                                           directed=FALSE)), is_true())
  g234 <- graph.atlas(234)
  expect_that(graph.isomorphic(g234, graph(c(1,6,2,6,3,6,4,6,5,6), n=7,
                                           directed=FALSE)), is_true())
})
