
context("graph.density")

test_that("graph.density works", {

  library(igraph)

  g <- erdos.renyi.game(50, 4/50)
  gd <- graph.density(g)
  gd2 <- ecount(g) / vcount(g) / (vcount(g)-1) * 2
  expect_that(gd, equals(gd2))

####

  g <- erdos.renyi.game(50, 4/50, dir=TRUE)
  gd <- graph.density(g)
  gd2 <- ecount(g) / vcount(g) / (vcount(g)-1)
  expect_that(gd, equals(gd2))

})
