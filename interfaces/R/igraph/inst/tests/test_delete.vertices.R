
context("delete.vertices")

test_that("delete.vertices works", {
  library(igraph)
  g <- graph.formula(A:B:C - D:E:F, D-E-F)
  
  g2 <- delete.vertices(g, "A")
  g3 <- delete.vertices(g, match("A", V(g)$name))
  
  expect_that(graph.isomorphic(g2, g3), is_true())
})
