
context("Various samplers")

test_that("Sampling from a Dirichlet works", {

  library(igraph)
  set.seed(42)
  sd <- sample_dirichlet(100, alpha=c(1, 1, 1))
  expect_that(dim(sd), equals(c(3, 100)))
  expect_that(colSums(sd), equals(rep(1, 100)))
  expect_that(mean(sd), equals(1/3))
  expect_that(sd(sd), equals(0.248901845755354))

  ## Corner cases
  sd1 <- sample_dirichlet(1, alpha=c(2, 2, 2))
  expect_that(dim(sd1), equals(c(3, 1)))
  sd0 <- sample_dirichlet(0, alpha=c(3, 3, 3))
  expect_that(dim(sd0), equals(c(3, 0)))

  ## Errors
  expect_that(sample_dirichlet(-1, alpha=c(1,1,1,1)),
              throws_error("should be non-negative"))
  expect_that(sample_dirichlet(5, alpha=c(1)),
              throws_error("must have at least two entries"))
  expect_that(sample_dirichlet(5, alpha=c(0, 1, 1)),
              throws_error("must be positive"))
  expect_that(sample_dirichlet(5, alpha=c(1, -1, -1)),
              throws_error("must be positive"))
})
