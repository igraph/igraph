
context("Fruchterman-Reingold layout")

test_that("", {

  library(igraph)
  set.seed(42)
  g <- graph.ring(10)
  l <- layout.fruchterman.reingold(g, niter=50)
  expect_that(sum(l), equals(10.7943032688805))

})
