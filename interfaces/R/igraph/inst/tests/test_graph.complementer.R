
context("complementer")

test_that("complementer works", {

  library(igraph)

  g <- sample_gnp(50, 3/50)
  g2 <- complementer(g)
  g3 <- complementer(g2)
  expect_that(graph.isomorphic(g, g3), is_true())
  
})


