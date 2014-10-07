
context("graph_from_adj_list")

test_that("graph_from_adj_list works", {

  library(igraph)
  
  g <- sample_gnp(100, 3/100)
  al <- as_adj_list(g)
  g2 <- graph_from_adj_list(al, mode="all")
  expect_that(graph.isomorphic(g, g2), is_true())

  ##
  
  g <- sample_gnp(100, 3/100, dir=TRUE)
  al <- as_adj_list(g, mode="out")
  g2 <- graph_from_adj_list(al, mode="out")
  expect_that(graph.isomorphic(g, g2), is_true())
})
