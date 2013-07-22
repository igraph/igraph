
context("Bonacich's power centrality")

test_that("Power centrality works", {
  library(igraph)

  ## Generate some test data from Bonacich, 1987:
  fig1 <- graph.formula( A -+ B -+ C:D )
  fig1.bp <- lapply(seq(0, 0.8, by=0.2), function(x)
                    round(bonpow(fig1, exponent=x), 2))
  expect_that(fig1.bp, equals(list(c(A=0.89, B=1.79, C=0, D=0),
                                   c(A=1.15, B=1.64, C=0, D=0),
                                   c(A=1.34, B=1.49, C=0, D=0),
                                   c(A=1.48, B=1.35, C=0, D=0),
                                   c(A=1.59, B=1.22, C=0, D=0))))

  g.c <- graph( c(1,2,1,3,2,4,3,5), dir=FALSE)
  bp.c <- lapply(seq(-.5, .5, by=0.1), function(x)
                 round(bonpow(g.c, exponent=x), 2)[c(1,2,4)])

  expect_that(bp.c, equals(list(c(0.00, 1.58, 0.00), c(0.73, 1.45, 0.36),
                                c(0.97, 1.34, 0.49), c(1.09, 1.27, 0.54),
                                c(1.15, 1.23, 0.58), c(1.20, 1.20, 0.60),
                                c(1.22, 1.17, 0.61), c(1.25, 1.16, 0.62),
                                c(1.26, 1.14, 0.63), c(1.27, 1.13, 0.64),
                                c(1.28, 1.12, 0.64))))

  g.d <- graph( c(1,2,1,3,1,4,2,5,3,6,4,7), dir=FALSE)
  bp.d <- lapply(seq(-.4, .4, by=0.1), function(x)
                 round(bonpow(g.d, exponent=x), 2)[c(1,2,5)])
  expect_that(bp.d, equals(list(c(1.62, 1.08, 0.54), c(1.62, 1.08, 0.54),
                                c(1.62, 1.08, 0.54), c(1.62, 1.08, 0.54),
                                c(1.62, 1.08, 0.54), c(1.62, 1.08, 0.54),
                                c(1.62, 1.08, 0.54), c(1.62, 1.08, 0.54),
                                c(1.62, 1.08, 0.54))))

  g.e <- graph( c(1,2,1,3,1,4,2,5,2,6,3,7,3,8,4,9,4,10), dir=FALSE)
  bp.e <- lapply(seq(-.4, .4, by=0.1), function(x)
                 round(bonpow(g.e, exponent=x), 2)[c(1,2,5)])
  expect_that(bp.e, equals(list(c(-1.00, 1.67, -0.33), c(0.36, 1.81, 0.12),
                                c( 1.00, 1.67,  0.33), c(1.30, 1.55, 0.43),
                                c( 1.46, 1.46,  0.49), c(1.57, 1.40, 0.52),
                                c( 1.63, 1.36,  0.54), c(1.68, 1.33, 0.56),
                                c( 1.72, 1.30,  0.57))))
  
  g.f <- graph( c(1,2,1,3,1,4,2,5,2,6,2,7,3,8,3,9,3,10,4,11,4,12,4,13),
               dir=FALSE)
  bp.f <- lapply(seq(-.4, .4, by=0.1), function(x)
                 round(bonpow(g.f, exponent=x), 2)[c(1,2,5)])
  expect_that(bp.f,
              equals(list(c(-1.72, 1.53, -0.57), c(-0.55, 2.03, -0.18),
                          c( 0.44, 2.05,  0.15), c( 1.01, 1.91,  0.34),
                          c( 1.33, 1.78,  0.44), c( 1.52, 1.67,  0.51),
                          c( 1.65, 1.59,  0.55), c( 1.74, 1.53,  0.58),
                          c( 1.80, 1.48,  0.60))))
})
