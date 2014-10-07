
context("edgelist")

test_that("as_edgelist works", {
  library(igraph)
  g <- sample_gnp(100, 3/100)
  e <- as_edgelist(g)
  g2 <- graph(t(e), n=vcount(g), dir=FALSE)
  expect_that(graph.isomorphic(g, g2), is_true())
})

