
context("layout.merge")

test_that("layout.merge works", {

  library(igraph)
  set.seed(42)

  g <- list(graph.ring(10), graph.ring(5))
  l <- lapply(g, layout.mds)
  l

  lm <- layout.merge(g, l)
  expect_that(is.matrix(lm), is_true())
  expect_that(ncol(lm), equals(2))
  expect_that(nrow(lm), equals(sum(sapply(g, vcount))))

##########

  ## Stress test
  for (i in 1:10) {
    g <- erdos.renyi.game(100, 2/100)
    l <- layout.mds(g)
    expect_that(dim(l), equals(c(vcount(g), 2)))
  }

})
