
context("girth")

test_that("girth works", {

  library(igraph)

  ## No circle in a tree
  g <- make_tree(1000, 3)
  gi <- girth(g)
  expect_that(gi$girth, equals(0))
  expect_that(as.vector(gi$circle), equals(numeric()))

  ## The worst case running time is for a ring
  g <- make_ring(100)
  gi <- girth(g)
  expect_that(gi$girth, equals(100))
  expect_that(sort(diff(as.vector(gi$circle))), equals(c(-99, rep(1, 98))))
})

