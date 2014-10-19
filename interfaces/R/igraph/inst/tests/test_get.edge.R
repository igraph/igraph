
context("ends")

test_that("ends works", {
  library(igraph)
  g <- sample_gnp(100, 3/100)
  edges <- unlist(lapply(seq_len(ecount(g)), ends, graph=g))
  g2 <- graph(edges, dir=FALSE, n=vcount(g))
  expect_that(graph.isomorphic(g, g2), is_true())
})
