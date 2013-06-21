
context("get.adjacency")

test_that("get.adjacency works", {

  library(igraph)

  g <- erdos.renyi.game(50, 1/50)
  A <- get.adjacency(g, sparse=FALSE)
  g2 <- graph.adjacency(A, mode="undirected")
  expect_that(graph.isomorphic(g, g2), is_true())

###

  A <- get.adjacency(g, sparse=TRUE)
  g2 <- graph.adjacency(A, mode="undirected")
  expect_that(graph.isomorphic(g, g2), is_true())

###

  g <- erdos.renyi.game(50, 2/50, directed=TRUE)
  A <- get.adjacency(g, sparse=FALSE)
  g2 <- graph.adjacency(A)
  expect_that(graph.isomorphic(g, g2), is_true())

###

  A <- get.adjacency(g, sparse=TRUE)
  g2 <- graph.adjacency(A)
  expect_that(graph.isomorphic(g, g2), is_true())

})
