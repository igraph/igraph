
context("graph.subisomorphic.vf2")

test_that("graph.subisomorphic.vf2 works", {

  library(igraph)
  set.seed(42)

  g1 <- erdos.renyi.game(20,6/20)
  g2 <- erdos.renyi.game(20,6/20)
  g <- g1 %du% g2

  ig1 <- graph.subisomorphic.vf2(g, g1)
  ig2 <- graph.subisomorphic.vf2(g, g2)

  expect_that(ig1$iso, is_true())
  expect_that(ig1$map12, equals(c(1:vcount(g1), rep(0, vcount(g2)))))
  expect_that(ig1$map21, equals(1:vcount(g1)))

  expect_that(ig2$iso, is_true())
  expect_that(ig2$map12, equals(c(rep(0, vcount(g1)), 1:vcount(g2))))
  expect_that(ig2$map21, equals(1:vcount(g2) + vcount(g1)))
  
})
