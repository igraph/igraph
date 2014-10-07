
context("compose")

test_that("compose works", {

  library(igraph)

  g1 <- sample_gnp(50, 3/50, directed=TRUE)
  gi <- graph( rep(1:vcount(g1), each=2), dir=TRUE )

  g2 <- compose(g1, gi)
  g3 <- compose(gi, g1)

  expect_that(graph.isomorphic(g1, g2), is_true())
  expect_that(graph.isomorphic(g1, g3), is_true())

})
