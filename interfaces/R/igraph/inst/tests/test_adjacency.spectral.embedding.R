
context("adjacency spectral embedding")

std <- function(x) {
  x <- zapsmall(x)
  apply(x, 2, function(col) {
    if (any(col < 0) && col[which(col != 0)[1]] < 0) { -col } else { col }
  })
}

mag_order <- function(x) {
  order(abs(x), sign(x), decreasing=TRUE)
}

mag_sort <- function(x) {
  x[mag_order(x)]
}

test_that("Undirected, unweighted case works", {
  library(igraph)
  set.seed(42)
  g <- random.graph.game(10, 15, type="gnm", directed=FALSE)

  no <- 7
  A <- g[]
  A <- A + 1/2 * diag(degree(g))
  ss <- eigen(A)

  U <- std(ss$vectors)
  X <- std(ss$vectors %*% sqrt(diag(abs(ss$values))))

  au_la <- embed_adjacency_matrix(g, no=no, which="la",
                                        cvec=degree(g)/2, scaled=TRUE)
  as_la <- embed_adjacency_matrix(g, no=no, which="la",
                                        cvec=degree(g)/2, scaled=FALSE)

  expect_that(as_la$D, equals(ss$values[1:no]))
  expect_that(au_la$D, equals(ss$values[1:no]))
  expect_that(std(as_la$X), equals(std(U[,1:no])))
  expect_that(std(au_la$X), equals(X[,1:no]))

  au_lm <- embed_adjacency_matrix(g, no=no, which="lm",
                                        cvec=degree(g)/2, scaled=TRUE)
  as_lm <- embed_adjacency_matrix(g, no=no, which="lm",
                                        cvec=degree(g)/2, scaled=FALSE)

  expect_that(as_lm$D, equals(mag_sort(ss$values)[1:no]))
  expect_that(au_lm$D, equals(mag_sort(ss$values)[1:no]))
  expect_that(std(as_lm$X), equals(std(U[,mag_order(ss$values)][,1:no])))
  expect_that(std(au_lm$X), equals(X[,mag_order(ss$values)][,1:no]))

  au_sa <- embed_adjacency_matrix(g, no=no, which="sa",
                                        cvec=degree(g)/2, scaled=TRUE)
  as_sa <- embed_adjacency_matrix(g, no=no, which="sa",
                                        cvec=degree(g)/2, scaled=FALSE)

  expect_that(as_sa$D, equals(ss$values[vcount(g)-1:no+1]))
  expect_that(au_sa$D, equals(ss$values[vcount(g)-1:no+1]))
  expect_that(std(as_sa$X), equals(std(U[,vcount(g)-1:no+1])))
  expect_that(std(au_sa$X), equals(X[,vcount(g)-1:no+1]))
})

test_that("Undirected, weighted case works", {
  library(igraph)
  set.seed(42)
  g <- random.graph.game(10, 20, type="gnm", directed=FALSE)
  E(g)$weight <- sample(1:5, ecount(g), replace=TRUE)

  no <- 3
  A <- g[]
  A <- A + 1/2 * diag(degree(g))
  ss <- eigen(A)

  U <- std(ss$vectors)
  X <- std(ss$vectors %*% sqrt(diag(abs(ss$values))))

  au_la <- embed_adjacency_matrix(g, no=no, which="la",
                                        cvec=degree(g)/2, scaled=TRUE)
  as_la <- embed_adjacency_matrix(g, no=no, which="la",
                                        cvec=degree(g)/2, scaled=FALSE)

  expect_that(as_la$D, equals(ss$values[1:no]))
  expect_that(std(as_la$X), equals(std(U[,1:no])))
  expect_that(au_la$D, equals(ss$values[1:no]))
  expect_that(std(au_la$X), equals(X[,1:no]))
  
  au_lm <- embed_adjacency_matrix(g, no=no, which="lm",
                                        cvec=degree(g)/2, scaled=TRUE)
  as_lm <- embed_adjacency_matrix(g, no=no, which="lm",
                                        cvec=degree(g)/2, scaled=FALSE)

  expect_that(as_lm$D, equals(mag_sort(ss$values)[1:no]))
  expect_that(au_lm$D, equals(mag_sort(ss$values)[1:no]))
  expect_that(std(as_lm$X), equals(std(U[,mag_order(ss$values)][,1:no])))
  expect_that(std(au_lm$X), equals(X[,mag_order(ss$values)][,1:no]))

  au_sa <- embed_adjacency_matrix(g, no=no, which="sa",
                                        cvec=degree(g)/2, scaled=TRUE)
  as_sa <- embed_adjacency_matrix(g, no=no, which="sa",
                                        cvec=degree(g)/2, scaled=FALSE)

  expect_that(std(as_sa$X), equals(std(U[,vcount(g)-1:no+1])))
  expect_that(std(au_sa$X), equals(X[,vcount(g)-1:no+1]))
})

test_that("Directed, unweighted case works", {
  library(igraph)
  set.seed(42)
  g <- random.graph.game(10, 20, type="gnm", directed=TRUE)

  no <- 3
  A <- g[]
  A <- A + 1/2 * diag(degree(g))
  ss <- svd(A)

  U <- std(ss$u)
  V <- std(ss$v)
  X <- std(ss$u %*% sqrt(diag(ss$d)))
  Y <- std(ss$v %*% sqrt(diag(ss$d)))

  au_la <- embed_adjacency_matrix(g, no=no, which="la",
                                        cvec=degree(g)/2, scaled=TRUE)
  as_la <- embed_adjacency_matrix(g, no=no, which="la",
                                        cvec=degree(g)/2, scaled=FALSE)

  expect_that(as_la$D, equals(ss$d[1:no]))
  expect_that(au_la$D, equals(ss$d[1:no]))
  expect_that(std(as_la$X), equals(std(U[,1:no])))
  expect_that(std(as_la$Y), equals(std(V[,1:no])))
  expect_that(std(au_la$X), equals(X[,1:no]))
  expect_that(std(au_la$Y), equals(Y[,1:no]))

  au_lm <- embed_adjacency_matrix(g, no=no, which="lm",
                                        cvec=degree(g)/2, scaled=TRUE)
  as_lm <- embed_adjacency_matrix(g, no=no, which="lm",
                                        cvec=degree(g)/2, scaled=FALSE)

  expect_that(as_lm$D, equals(ss$d[1:no]))
  expect_that(au_lm$D, equals(ss$d[1:no]))
  expect_that(std(as_lm$X), equals(std(U[,1:no])))
  expect_that(std(as_lm$Y), equals(std(V[,1:no])))
  expect_that(std(au_lm$X), equals(X[,1:no]))
  expect_that(std(au_lm$Y), equals(Y[,1:no]))

  au_sa <- embed_adjacency_matrix(g, no=no, which="sa",
                                        cvec=degree(g)/2, scaled=TRUE)
  as_sa <- embed_adjacency_matrix(g, no=no, which="sa",
                                        cvec=degree(g)/2, scaled=FALSE)

  expect_that(as_sa$D, equals(ss$d[vcount(g)-1:no+1]))
  expect_that(au_sa$D, equals(ss$d[vcount(g)-1:no+1]))
  expect_that(std(as_sa$X), equals(std(U[,vcount(g)-1:no+1])))
  expect_that(std(as_sa$Y), equals(std(V[,vcount(g)-1:no+1])))
  expect_that(std(au_sa$X), equals(X[,vcount(g)-1:no+1]))
  expect_that(std(au_sa$Y), equals(Y[,vcount(g)-1:no+1]))
})

test_that("Directed, weighted case works", {
  library(igraph)
  set.seed(42)
  g <- random.graph.game(10, 20, type="gnm", directed=TRUE)
  E(g)$weight <- sample(1:5, ecount(g), replace=TRUE)

  no <- 3
  A <- g[]
  A <- A + 1/2 * diag(degree(g))
  ss <- svd(A)

  U <- std(ss$u)
  V <- std(ss$v)
  X <- std(ss$u %*% sqrt(diag(ss$d)))
  Y <- std(ss$v %*% sqrt(diag(ss$d)))

  au_la <- embed_adjacency_matrix(g, no=no, which="la",
                                        cvec=degree(g)/2, scaled=TRUE)
  as_la <- embed_adjacency_matrix(g, no=no, which="la",
                                        cvec=degree(g)/2, scaled=FALSE)

  expect_that(std(as_la$X), equals(std(U[,1:no])))
  expect_that(std(as_la$Y), equals(std(V[,1:no])))
  expect_that(std(au_la$X), equals(X[,1:no]))
  expect_that(std(au_la$Y), equals(Y[,1:no]))

  au_lm <- embed_adjacency_matrix(g, no=no, which="lm",
                                        cvec=degree(g)/2, scaled=TRUE)
  as_lm <- embed_adjacency_matrix(g, no=no, which="lm",
                                        cvec=degree(g)/2, scaled=FALSE)

  expect_that(std(as_lm$X), equals(std(U[,1:no])))
  expect_that(std(as_lm$Y), equals(std(V[,1:no])))
  expect_that(std(au_lm$X), equals(X[,1:no]))
  expect_that(std(au_lm$Y), equals(Y[,1:no]))

  au_sa <- embed_adjacency_matrix(g, no=no, which="sa",
                                        cvec=degree(g)/2, scaled=TRUE)
  as_sa <- embed_adjacency_matrix(g, no=no, which="sa",
                                        cvec=degree(g)/2, scaled=FALSE)

  expect_that(std(as_sa$X), equals(std(U[,vcount(g)-1:no+1])))
  expect_that(std(as_sa$Y), equals(std(V[,vcount(g)-1:no+1])))
  expect_that(std(au_sa$X), equals(X[,vcount(g)-1:no+1]))
  expect_that(std(au_sa$Y), equals(Y[,vcount(g)-1:no+1]))
})

test_that("Issue #50 is resolved", {

  library(igraph)
  set.seed(12345)

  g <- erdos.renyi.game(15, .4)
  w <- -log(runif(ecount(g)))
  X1 <- embed_adjacency_matrix(g, 2, weights= w)

  E(g)$weight <- w
  X2 <- embed_adjacency_matrix(g, 2)

  expect_that(X1$D, equals(X2$D))

})

test_that("Issue #51 is resolved", {
  library(igraph)
  set.seed(12345)

  pref.matrix <- diag(0.2, 2) + 0.2
  block.sizes <- c(800, 800)
  n <- sum(block.sizes)
  g <- sbm.game(n, pref.matrix, block.sizes, directed=TRUE)

  for (i in 1:25) {
    ase <-  embed_adjacency_matrix(g, 2)
    expect_that(mean(ase$X %*% t(ase$Y)), equals(0.299981018354173))
  }
})
