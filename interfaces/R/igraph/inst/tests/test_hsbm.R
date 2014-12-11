
context("Hierarchical stochastic block models")

test_that("HSBM works", {
  library(igraph)
  set.seed(42)

  C <- matrix(c(1  , 1/2,   0,
                1/2,   0, 1/2,
                0  , 1/2, 1/2), nrow=3)
                
  g <- sample_hierarchical_sbm(100, 10, rho=c(3,3,4)/10, C=C, p=0)
  expect_that(ecount(g), equals(172))
  expect_that(vcount(g), equals(100))
  expect_that(is.directed(g), is_false())

  set.seed(42)

  g2 <- sample_hierarchical_sbm(100, 10, rho=c(3,3,4)/10, C=C, p=1)
  expect_that(ecount(g2), equals(ecount(g) + 10 * 9 * (90 + 10) / 2))
  expect_that(vcount(g2), equals(100))
  expect_that(is.simple(g2), is_true())

  set.seed(42)
  
  g3 <- sample_hierarchical_sbm(100, 10, rho=c(3,3,4)/10, C=C, p=1e-15)
  expect_that(ecount(g3), equals(ecount(g)))
  expect_that(vcount(g3), equals(100))
  expect_that(is.simple(g3), is_true())

  set.seed(42)
  
  g4 <- sample_hierarchical_sbm(100, 10, rho=c(3,3,4)/10, C=C, p=1-1e-15)
  expect_that(ecount(g4), equals(ecount(g2)))
  expect_that(vcount(g4), equals(100))
  expect_that(is.simple(g4), is_true())

})

test_that("HSBM with 1 cluster per block works", {
  library(igraph)
  res <- Matrix(0, nrow=10, ncol=10)
  res[6:10, 1:5] <- res[1:5, 6:10] <- 1
  g <- sample_hierarchical_sbm(10, 5, rho=1, C=matrix(0), p=1)
  expect_that(g[], equals(res))
})

test_that("HSBM with list arguments works", {
  library(igraph)

  b <- 5
  C <- matrix(c(1  , 1/2,   0,
                1/2,   0, 1/2,
                0  , 1/2, 1/2), nrow=3)
  m <- 10
  rho <- c(3,3,4)/10

  set.seed(42)
  g <- sample_hierarchical_sbm(b*m, m, rho=rho, C=C, p=0)

  set.seed(42)
  g2 <- sample_hierarchical_sbm(b*m, rep(m, b), rho=rho, C=C, p=0)
  expect_that(g[], equals(g2[]))

  set.seed(42)
  g3 <- sample_hierarchical_sbm(b*m, m, rho=replicate(b, rho, simplify=FALSE), C=C, p=0)
  expect_that(g[], equals(g3[]))

  set.seed(42)
  g4 <- sample_hierarchical_sbm(b*m, m, rho=rho, C=replicate(b, C, simplify=FALSE), p=0)
  expect_that(g[], equals(g4[]))

  expect_that(sample_hierarchical_sbm(b*m, rep(m, b), rho=list(rho, rho), C=C, p=0),
              throws_error("Lengths of `m', `rho' and `C' must match"))

  ###

  n <- function(x) x/sum(x)
  
  rho1 <- n(c(1,2))
  C1 <- matrix(0, nrow=2, ncol=2)
  rho2 <- n(c(3,3,4))
  C2 <- matrix(0, nrow=3, ncol=3)
  rho3 <- 1
  C3 <- matrix(0)
  rho4 <- n(c(2,1))
  C4 <- matrix(0, nrow=2, ncol=2)

  gg1 <- sample_hierarchical_sbm(21, m=c(3, 10, 5, 3), rho=list(rho1, rho2, rho3, rho4),
                   C=list(C1, C2, C3, C4), p=1)
  expect_that(is.simple(gg1), is_true())

  set.seed(42)
  gg11 <- sample_hierarchical_sbm(21, m=c(3, 10, 5, 3), rho=list(rho1, rho2, rho3, rho4),
                    C=list(C1, C2, C3, C4), p=1-1e-10)
  expect_that(gg1[], equals(gg11[]))

  rho1 <- n(c(1,2))
  C1 <- matrix(1, nrow=2, ncol=2)
  rho2 <- n(c(3,3,4))
  C2 <- matrix(1, nrow=3, ncol=3)
  rho3 <- 1
  C3 <- matrix(1)
  rho4 <- n(c(2,1))
  C4 <- matrix(1, nrow=2, ncol=2)
  gg2 <- sample_hierarchical_sbm(21, m=c(3, 10, 5, 3), rho=list(rho1, rho2, rho3, rho4),
                   C=list(C1, C2, C3, C4), p=0)
  expect_that(is.simple(gg2), is_true())

  gg22 <- sample_hierarchical_sbm(21, m=c(3, 10, 5, 3), rho=list(rho1, rho2, rho3, rho4),
                   C=list(C1, C2, C3, C4), p=1)
  expect_that(gg1[] + gg2[], equals(gg22[]))
})
