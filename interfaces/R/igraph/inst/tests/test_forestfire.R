
context("forest.fire.game")

test_that("forest.fire.game works", {
  
  library(igraph)
  set.seed(42)

  pars <- list(sparse=c(0.35, 0.2/0.35),
               densifying=c(0.37, 0.32/0.37),
               dense=c(0.38, 0.38/0.37))

  N <- 5000
  G <- lapply(pars, function(x) forest.fire.game(N, fw=x[1], bw=x[2]))
  
  xv <- log(2:N)
  
  co <- sapply(G, function(x) {
    yv <- log(cumsum(degree(x, mode="out"))[-1])
    coef(lm( yv ~ xv ))[2]
  })
  
  expect_that(co, equals(structure(c(1.06045500245466,
                                     1.22800967143684,
                                     1.96234121488344),
                                   .Names = c("sparse.xv",
                                     "densifying.xv", "dense.xv"))))

})
