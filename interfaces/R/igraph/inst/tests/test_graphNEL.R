
context("graphNEL conversion")

test_that("graphNEL conversion works", {

  library(igraph)
  library(graph, warn.conflicts=FALSE)

  g <- erdos.renyi.game(100, 5/100)
  N <- igraph.to.graphNEL(g)
  g2 <- igraph.from.graphNEL(N)
  gi <- graph.isomorphic.vf2(g, g2)
  expect_that(gi$iso, is_true())
  expect_that(gi$map12, equals(1:vcount(g)))
  expect_that(gi$map21, equals(1:vcount(g)))

  ## Attributes

  V(g)$name <- as.character(vcount(g):1)
  E(g)$weight <- sample(1:10, ecount(g), replace=TRUE)
  g$name <- "Foobar"

  N <- igraph.to.graphNEL(g)
  g2 <- igraph.from.graphNEL(N)
  expect_that(graph.isomorphic(g, g2), is_true())
  expect_that(V(g)$name, equals(V(g2)$name))

  A <- get.adjacency(g, attr="weight", sparse=FALSE)
  A2 <- get.adjacency(g2, attr="weight", sparse=FALSE)
  expect_that(A, equals(A))
  expect_that(g$name, equals(g2$name))

  suppressWarnings(unloadNamespace("graph"))
})
