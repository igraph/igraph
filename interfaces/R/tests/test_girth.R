
context("girth")

test_that("girth works", {

  library(igraph)

  ## No circle in a tree
  g <- graph.tree(1000, 3)
  gi <- girth(g)
  expect_that(gi$girth, equals(0))
  expect_that(gi$circle, equals(numeric()))

  ## The worst case running time is for a ring
  g <- graph.ring(100)
  gi <- girth(g)
  expect_that(gi$girth, equals(100))
  expect_that(sort(diff(gi$circle)), equals(c(-99, rep(1, 98))))
})

