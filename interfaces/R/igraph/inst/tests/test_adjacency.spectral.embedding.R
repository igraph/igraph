
context("adjacency spectral embedding")

test_that("adjacency.spectral.embedding works", {
  library(igraph)

  c <- 1/2
  no <- 3
  
  g <- graph.tree(10, 3, mode="out")
  asm <- adjacency.spectral.embedding(g, no=no, c=c)

  A <- get.adjacency(g)
  A <- A + c * diag(degree(g))
  asm2 <- svd(A)

  std <- function(x) {
    x <- zapsmall(x)
    apply(x, 2, function(col) {
      if (any(col < 0) && col[which(col != 0)[1]] < 0) { -col } else { col }
    })
  }
  
  expect_that(asm$D, equals(asm2$d[1:no]))
  expect_that(std(asm$U), equals(std(asm2$u[,1:no])))
  expect_that(std(asm$V), equals(std(asm2$v[,1:no])))
  
})
