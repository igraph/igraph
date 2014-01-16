
context("SIR epidemics model on a network")

test_that("SIR works", {

  set.seed(42)
  library(digest)
  library(igraph)

  g <- erdos.renyi.game(50, 50, type="gnm")
  res <- sir(g, beta=5, gamma=1, no.sim=10)
  res2 <- list()
  for (i in 1:10) {
    res2[[i]] <- list(times=res$times[[i]], NS=res$no_s[[i]],
                      NI=res$no_i[[i]], NR=res$no_r[[i]])
  }
  class(res2) <- "sir"
  expect_that(digest(res2), equals("bc42d0cbe0bb3321e83979c0432f9cea"))
})

