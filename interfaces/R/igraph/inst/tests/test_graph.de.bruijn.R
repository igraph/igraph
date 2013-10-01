
context("graph.de.bruijn")

test_that("graph.de.bruijn works", {

  library(igraph)
  g <- graph.de.bruijn(2,1)
  g2 <- graph.de.bruijn(2,2)
  g3 <- line.graph(g)

  expect_that(graph.isomorphic(g3, graph(c(1,1,3,1,1,2,3,2,2,3,
                                           4,3,2,4,4,4))), is_true())
  expect_that(graph.isomorphic(g2, g3), is_true())
  
})
