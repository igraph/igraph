
context("decompose.graph")

test_that("decompose.graph works", {
  library(igraph)
  g <- erdos.renyi.game(1000, 1/1500)
  G <- decompose.graph(g)
  clu <- clusters(g)
  Gsizes <- sapply(G, vcount)
  expect_that(sort(clu$csize), equals(sort(Gsizes)))
})
