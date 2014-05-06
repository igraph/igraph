
context("adjacency spectral embedding")

test_that("adjacency.spectral.embedding works", {
  library(igraph)

  no <- 3
  g <- graph.tree(10, 3, mode="out")
  asm <- adjacency.spectral.embedding(g, no=no, which="la", cvec=degree(g)/2)
  asmu <- adjacency.spectral.embedding(g, no=no, which="la", cvec=degree(g)/2,
                                       scaled=FALSE)

  A <- get.adjacency(g)
  A <- A + 1/2 * diag(degree(g))
  asm2 <- svd(A)

  std <- function(x) {
    x <- zapsmall(x)
    apply(x, 2, function(col) {
      if (any(col < 0) && col[which(col != 0)[1]] < 0) { -col } else { col }
    })
  }

  expect_that(std(asmu$X), equals(std(asm2$u[,1:no])))
  expect_that(std(asmu$Y), equals(std(asm2$v[,1:no])))
  
  X <- asm2$u[,1:no] %*% sqrt(diag(asm2$d[1:no]))
  Y <- asm2$v[,1:no] %*% sqrt(diag(asm2$d[1:no]))
  
  expect_that(std(asm$X), equals(std(X)))
  expect_that(std(asm$Y), equals(std(Y)))

})

test_that("weighted adjacency spectral embedding works", {
  library(igraph)

  set.seed(42)
  no <- 3 
  g <- graph.tree(10, 3, mode="out")
  E(g)$weight <- sample(1:5, ecount(g), replace=TRUE)
  asm <- adjacency.spectral.embedding(g, no=no, which="la", cvec=degree(g)/2)
  asmu <- adjacency.spectral.embedding(g, no=no, which="la", cvec=degree(g)/2,
                                       scaled=FALSE)

  A <- g[]
  A <- A + 1/2 * diag(degree(g))
  asm2 <- svd(A)

  std <- function(x) {
    x <- zapsmall(x)
    apply(x, 2, function(col) {
      if (any(col < 0) && col[which(col != 0)[1]] < 0) { -col } else { col }
    })
  }

  expect_that(std(asmu$X), equals(std(asm2$u[,1:no])))
  expect_that(std(asmu$Y), equals(std(asm2$v[,1:no])))
  
  X <- asm2$u[,1:no] %*% sqrt(diag(asm2$d[1:no]))
  Y <- asm2$v[,1:no] %*% sqrt(diag(asm2$d[1:no]))
  
  expect_that(std(asm$X), equals(std(X)))
  expect_that(std(asm$Y), equals(std(Y)))

})
