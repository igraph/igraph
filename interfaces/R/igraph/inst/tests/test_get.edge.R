
context("get.edge")

test_that("get.edge works", {
  library(igraph)
  g <- erdos.renyi.game(100, 3/100)
  edges <- unlist(lapply(seq_len(ecount(g)), get.edge, graph=g))
  g2 <- graph(edges, dir=FALSE, n=vcount(g))
  expect_that(graph.isomorphic(g, g2), is_true())
})
