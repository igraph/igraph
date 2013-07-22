
context("watts.strogatz.game")

test_that("watts.strogatz.game works", {

  library(igraph)

  for (i in 1:50) {
    p <- runif(1)
    d <- sample(1:3, 1)
    nei <- sample(2:5, 1)
    g <- watts.strogatz.game(d, 10, nei, p, loops=FALSE)
    expect_that(any(is.loop(g)), is_false())
  }

})
