
context("Eigenproblems")

test_that("spectrum works for symmetric matrices", {
  library(igraph)
  set.seed(42)

  std <- function(x) {
    x <- zapsmall(x)
    apply(x, 2, function(col) {
      if (any(col < 0) && col[which(col != 0)[1]] < 0) { -col } else { col }
    })
  }
  
  g <- sample_gnp(50, 5/50)
  e0 <- eigen(as_adj(g, sparse=FALSE))

  e1 <- spectrum(g, which=list(howmany=4, pos="LA"))
  expect_that(e0$values[1:4], equals(e1$values))
  expect_that(std(e0$vectors[,1:4]), equals(std(e1$vectors)))

  e2 <- spectrum(g, which=list(howmany=4, pos="SA"))
  expect_that(e0$values[50:47], equals(e2$values))
  expect_that(std(e0$vectors[,50:47]), equals(std(e2$vectors)))

})
