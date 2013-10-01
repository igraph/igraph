
context("neighbors")

test_that("neighbors works", {

  library(igraph)

  g <- erdos.renyi.game(100, 20/100)
  al <- get.adjlist(g, mode="all")
  for (i in 1:length(al)) {
    n <- neighbors(g, i, mode="out")
    expect_that(sort(n), equals(al[[i]]))
  }

})
