
context("SIR epidemics model on a network")

test_that("SIR works", {

  set.seed(42)
  library(digest)
  library(igraph)

  g <- erdos.renyi.game(50, 50, type="gnm")
  res <- sir(g, beta=5, gamma=1, no.sim=10)
  expect_that(digest(res), equals("bc42d0cbe0bb3321e83979c0432f9cea"))
})

