
context("Hierarchical stochastic block models")

test_that("HSBM works", {
  library(igraph)
  set.seed(42)

  C <- matrix(c(1  , 1/2,   0,
                1/2,   0, 1/2,
                0  , 1/2, 1/2), nrow=3)
                
  g <- hsbm.game(100, 10, rho=c(3,3,4)/10, C=C, p=0)
  expect_that(ecount(g), equals(172))
  expect_that(vcount(g), equals(100))
  expect_that(is.directed(g), is_false())

  set.seed(42)

  g2 <- hsbm.game(100, 10, rho=c(3,3,4)/10, C=C, p=1)
  expect_that(ecount(g2), equals(ecount(g) + 10 * 9 * (90 + 10) / 2))
  expect_that(vcount(g2), equals(100))
  expect_that(is.simple(g2), is_true())

  set.seed(42)
  
  g3 <- hsbm.game(100, 10, rho=c(3,3,4)/10, C=C, p=1e-15)
  expect_that(ecount(g3), equals(ecount(g)))
  expect_that(vcount(g3), equals(100))
  expect_that(is.simple(g3), is_true())

  set.seed(42)
  
  g4 <- hsbm.game(100, 10, rho=c(3,3,4)/10, C=C, p=1-1e-15)
  expect_that(ecount(g4), equals(ecount(g2)))
  expect_that(vcount(g4), equals(100))
  expect_that(is.simple(g4), is_true())

})
