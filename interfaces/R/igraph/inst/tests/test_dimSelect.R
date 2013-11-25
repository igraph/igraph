
context("Dimensionality selection")

test_that("dimensionality selection works", {
  library(igraph)
  set.seed(42)

  k <- graph.famous("zachary")
  ev <- eigen(get.adjacency(k), only.values=TRUE)$values
  kdim <- dimSelect(ev)
  expect_that(kdim, equals(4))

  expect_that(dimSelect(1:100), equals(50))

  ## Some regression tests
  expect_that(dimSelect(runif(100)), equals(69))
  expect_that(dimSelect(runif(100)), equals(88))
  expect_that(dimSelect(runif(100)), equals(3))
  expect_that(dimSelect(runif(100)), equals(99))
  
})
