
context("graph.formula")

test_that("simplify argument works", {

  library(igraph)
  g1 <- graph.formula(1-1, 1-2, 1-2)
  g2 <- graph.formula(1-1, 1-2, 1-2, simplify=FALSE)

  expect_that(vcount(g1), equals(2))
  expect_that(ecount(g1), equals(1))

  expect_that(vcount(g2), equals(2))
  expect_that(ecount(g2), equals(3))
  
})
