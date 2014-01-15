
context("SIR epidemics model on a network")

test_that("SIR works", {

  set.seed(42)
  library(digest)
  library(igraph)

  g <- erdos.renyi.game(50, 50, type="gnm")
  res <- sir(g, beta=5, gamma=1, no_sim=10)
  expect_that(digest(res), equals("1fb77ab6f4184333470852028edc4f5e"))
})

