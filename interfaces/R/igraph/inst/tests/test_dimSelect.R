
context("Dimensionality selection")

test_that("dimensionality selection works", {
  library(igraph)
  set.seed(42)

  k <- graph.famous("zachary")
  ev <- eigen(get.adjacency(k), only.values=TRUE)$values
  kdim <- dim_select(ev)
  expect_that(kdim, equals(4))

  expect_that(dim_select(1:100), equals(50))

  ## Some regression tests
  expect_that(dim_select(runif(100)), equals(69))
  expect_that(dim_select(runif(100)), equals(88))
  expect_that(dim_select(runif(100)), equals(3))
  expect_that(dim_select(runif(100)), equals(99))

  ## Some more meaningful tests
  x <- c(rnorm(50, mean=0, sd=1), rnorm(50, mean=5, sd=1))
  expect_that(dim_select(x), equals(50))

  x <- c(rnorm(10, mean=0, sd=1), rnorm(90, mean=2, sd=1))
  expect_that(dim_select(x), equals(10))

  x <- c(10, rnorm(99, mean=0, sd=1))
  expect_that(dim_select(x), equals(1))

  x <- c(rnorm(99, mean=0, sd=1), 10)
  expect_that(dim_select(x), equals(99))
})
