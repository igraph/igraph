
context("Spectral embedding of the Laplacian")

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

test_that("Undirected, unweighted, D-A case works", {
  library(igraph)
  set.seed(42)
  g <- random.graph.game(10, 20, type="gnm", directed=FALSE)

  no <- 3
  A <- diag(degree(g)) - g[]
  ss <- eigen(A)

  D <- ss$values
  U <- ss$vectors
  X <- std(ss$vectors %*% sqrt(diag(ss$values)))
  Y <- std(ss$vectors %*% sqrt(diag(ss$values)))

  ## LA

  au_la <- laplacian.spectral.embedding(g, no=no, which="la", type="D-A",
                                        scaled=TRUE)
  as_la <- laplacian.spectral.embedding(g, no=no, which="la", type="D-A",
                                        scaled=FALSE)

  expect_that(au_la$D, equals(D[1:no]))
  expect_that(std(au_la$X), equals(std(X[,1:no])))
  expect_that(as_la$D, equals(D[1:no]))
  expect_that(std(as_la$X), equals(std(U[,1:no])))

  ## LM

  au_lm <- laplacian.spectral.embedding(g, no=no, which="lm", type="D-A",
                                        scaled=TRUE)
  as_lm <- laplacian.spectral.embedding(g, no=no, which="lm", type="D-A",
                                        scaled=FALSE)

  expect_that(au_lm$D, equals(mag_sort(D)[1:no]))
  expect_that(std(au_lm$X), equals(std(X[,mag_order(D)][,1:no])))
  expect_that(as_lm$D, equals(mag_sort(D)[1:no]))
  expect_that(std(as_lm$X), equals(std(U[,mag_order(D)][,1:no])))

  ## SA

  au_sa <- laplacian.spectral.embedding(g, no=no, which="sa", type="D-A",
                                        scaled=TRUE)
  as_sa <- laplacian.spectral.embedding(g, no=no, which="sa", type="D-A",
                                        scaled=FALSE)

  expect_that(au_sa$D, equals(D[vcount(g)-1:no+1]))
  expect_that(std(au_sa$X), equals(std(X[,vcount(g)-1:no+1])))
  expect_that(as_sa$D, equals(D[vcount(g)-1:no+1]))
  expect_that(std(as_sa$X), equals(std(U[,vcount(g)-1:no+1])))

})

test_that("Undirected, unweighted, DAD case works", {
  library(igraph)
  set.seed(42)
  g <- random.graph.game(10, 20, type="gnm", directed=FALSE)

  no <- 3
  D12 <- diag(1/sqrt(degree(g)))
  A <- D12 %*% g[] %*% D12
  ss <- eigen(A)

  D <- ss$values
  U <- ss$vectors
  X <- std(ss$vectors %*% sqrt(diag(abs(ss$values))))
  Y <- std(ss$vectors %*% sqrt(diag(abs(ss$values))))

  ## LA

  au_la <- laplacian.spectral.embedding(g, no=no, which="la", type="DAD",
                                        scaled=TRUE)
  as_la <- laplacian.spectral.embedding(g, no=no, which="la", type="DAD",
                                        scaled=FALSE)

  expect_that(au_la$D, equals(abs(D[1:no])))
  expect_that(std(au_la$X), equals(std(X[,1:no])))
  expect_that(as_la$D, equals(D[1:no]))
  expect_that(std(as_la$X), equals(std(U[,1:no])))

  ## LM

  au_lm <- laplacian.spectral.embedding(g, no=no, which="lm", type="DAD",
                                        scaled=TRUE)
  as_lm <- laplacian.spectral.embedding(g, no=no, which="lm", type="DAD",
                                        scaled=FALSE)

  expect_that(au_lm$D, equals(mag_sort(D)[1:no]))
  expect_that(std(au_lm$X), equals(std(X[,mag_order(D)][,1:no])))
  expect_that(as_lm$D, equals(mag_sort(D)[1:no]))
  expect_that(std(as_lm$X), equals(std(U[,mag_order(D)][,1:no])))

  ## SA

  au_sa <- laplacian.spectral.embedding(g, no=no, which="sa", type="DAD",
                                        scaled=TRUE)
  as_sa <- laplacian.spectral.embedding(g, no=no, which="sa", type="DAD",
                                        scaled=FALSE)

  expect_that(au_sa$D, equals(D[vcount(g)-1:no+1]))
  expect_that(std(au_sa$X), equals(std(X[,vcount(g)-1:no+1])))
  expect_that(as_sa$D, equals(D[vcount(g)-1:no+1]))
  expect_that(std(as_sa$X), equals(std(U[,vcount(g)-1:no+1])))
})

test_that("Undirected, unweighted, I-DAD case works", {
  library(igraph)
  set.seed(42)
  g <- random.graph.game(10, 20, type="gnm", directed=FALSE)

  no <- 3
  D12 <- diag(1/sqrt(degree(g)))
  A <- diag(vcount(g)) - D12 %*% g[] %*% D12
  ss <- eigen(A)

  D <- ss$values
  U <- ss$vectors
  X <- std(ss$vectors %*% sqrt(diag(abs(ss$values))))
  Y <- std(ss$vectors %*% sqrt(diag(abs(ss$values))))

  ## LA

  au_la <- laplacian.spectral.embedding(g, no=no, which="la", type="I-DAD",
                                        scaled=TRUE)
  as_la <- laplacian.spectral.embedding(g, no=no, which="la", type="I-DAD",
                                        scaled=FALSE)

  expect_that(au_la$D, equals(abs(D[1:no])))
  expect_that(std(au_la$X), equals(std(X[,1:no])))
  expect_that(as_la$D, equals(D[1:no]))
  expect_that(std(as_la$X), equals(std(U[,1:no])))

  ## LM

  au_lm <- laplacian.spectral.embedding(g, no=no, which="lm", type="I-DAD",
                                        scaled=TRUE)
  as_lm <- laplacian.spectral.embedding(g, no=no, which="lm", type="I-DAD",
                                        scaled=FALSE)

  expect_that(au_lm$D, equals(mag_sort(D)[1:no]))
  expect_that(std(au_lm$X), equals(std(X[,mag_order(D)][,1:no])))
  expect_that(as_lm$D, equals(mag_sort(D)[1:no]))
  expect_that(std(as_lm$X), equals(std(U[,mag_order(D)][,1:no])))

  ## SA

  au_sa <- laplacian.spectral.embedding(g, no=no, which="sa", type="I-DAD",
                                        scaled=TRUE)
  as_sa <- laplacian.spectral.embedding(g, no=no, which="sa", type="I-DAD",
                                        scaled=FALSE)

  expect_that(au_sa$D, equals(D[vcount(g)-1:no+1]))
  expect_that(std(au_sa$X), equals(std(X[,vcount(g)-1:no+1])))
  expect_that(as_sa$D, equals(D[vcount(g)-1:no+1]))
  expect_that(std(as_sa$X), equals(std(U[,vcount(g)-1:no+1])))
})

test_that("Undirected, weighted, D-A case works", {
  library(igraph)
  set.seed(42*42)
  g <- random.graph.game(10, 20, type="gnm", directed=FALSE)
  E(g)$weight <- sample(1:5, ecount(g), replace=TRUE)

  no <- 3
  A <- diag(graph.strength(g)) - g[]
  ss <- eigen(A)

  D <- ss$values
  U <- ss$vectors
  X <- std(ss$vectors %*% sqrt(diag(abs(D))))
  Y <- std(ss$vectors %*% sqrt(diag(abs(D))))

  ## LA

  au_la <- laplacian.spectral.embedding(g, no=no, which="la", type="D-A",
                                        scaled=TRUE)
  as_la <- laplacian.spectral.embedding(g, no=no, which="la", type="D-A",
                                        scaled=FALSE)

  expect_that(au_la$D, equals(abs(D[1:no])))
  expect_that(std(au_la$X), equals(std(X[,1:no])))
  expect_that(as_la$D, equals(D[1:no]))
  expect_that(std(as_la$X), equals(std(U[,1:no])))

  ## LM

  au_lm <- laplacian.spectral.embedding(g, no=no, which="lm", type="D-A",
                                        scaled=TRUE)
  as_lm <- laplacian.spectral.embedding(g, no=no, which="lm", type="D-A",
                                        scaled=FALSE)

  expect_that(au_lm$D, equals(mag_sort(D)[1:no]))
  expect_that(std(au_lm$X), equals(std(X[,mag_order(D)][,1:no])))
  expect_that(as_lm$D, equals(mag_sort(D)[1:no]))
  expect_that(std(as_lm$X), equals(std(U[,mag_order(D)][,1:no])))

  ## SA

  au_sa <- laplacian.spectral.embedding(g, no=no, which="sa", type="D-A",
                                        scaled=TRUE)
  as_sa <- laplacian.spectral.embedding(g, no=no, which="sa", type="D-A",
                                        scaled=FALSE)

  expect_that(au_sa$D, equals(D[vcount(g)-1:no+1]))
  expect_that(std(au_sa$X), equals(X[,vcount(g)-1:no+1],
                                   tolerance=.Machine$double.eps ^ 0.25))
  expect_that(as_sa$D, equals(D[vcount(g)-1:no+1]))
  expect_that(std(as_sa$X), equals(std(U[,vcount(g)-1:no+1])))

})

test_that("Undirected, unweighted, DAD case works", {
  library(igraph)
  set.seed(42)
  g <- random.graph.game(10, 20, type="gnm", directed=FALSE)

  no <- 3
  D12 <- diag(1/sqrt(degree(g)))
  A <- D12 %*% g[] %*% D12
  ss <- eigen(A)

  D <- ss$values
  U <- ss$vectors
  X <- std(ss$vectors %*% sqrt(diag(abs(ss$values))))
  Y <- std(ss$vectors %*% sqrt(diag(abs(ss$values))))

  ## LA

  au_la <- laplacian.spectral.embedding(g, no=no, which="la", type="DAD",
                                        scaled=TRUE)
  as_la <- laplacian.spectral.embedding(g, no=no, which="la", type="DAD",
                                        scaled=FALSE)

  expect_that(au_la$D, equals(abs(D[1:no])))
  expect_that(std(au_la$X), equals(std(X[,1:no])))
  expect_that(as_la$D, equals(D[1:no]))
  expect_that(std(as_la$X), equals(std(U[,1:no])))

  ## LM

  au_lm <- laplacian.spectral.embedding(g, no=no, which="lm", type="DAD",
                                        scaled=TRUE)
  as_lm <- laplacian.spectral.embedding(g, no=no, which="lm", type="DAD",
                                        scaled=FALSE)

  expect_that(au_lm$D, equals(mag_sort(D)[1:no]))
  expect_that(std(au_lm$X), equals(std(X[,mag_order(D)][,1:no])))
  expect_that(as_lm$D, equals(mag_sort(D)[1:no]))
  expect_that(std(as_lm$X), equals(std(U[,mag_order(D)][,1:no])))

  ## SA

  au_sa <- laplacian.spectral.embedding(g, no=no, which="sa", type="DAD",
                                        scaled=TRUE)
  as_sa <- laplacian.spectral.embedding(g, no=no, which="sa", type="DAD",
                                        scaled=FALSE)

  expect_that(au_sa$D, equals(D[vcount(g)-1:no+1]))
  expect_that(std(au_sa$X), equals(std(X[,vcount(g)-1:no+1])))
  expect_that(as_sa$D, equals(D[vcount(g)-1:no+1]))
  expect_that(std(as_sa$X), equals(std(U[,vcount(g)-1:no+1])))
})

test_that("Undirected, unweighted, I-DAD case works", {
  library(igraph)
  set.seed(42)
  g <- random.graph.game(10, 20, type="gnm", directed=FALSE)

  no <- 3
  D12 <- diag(1/sqrt(degree(g)))
  A <- diag(vcount(g)) - D12 %*% g[] %*% D12
  ss <- eigen(A)

  D <- ss$values
  U <- ss$vectors
  X <- std(ss$vectors %*% sqrt(diag(abs(ss$values))))
  Y <- std(ss$vectors %*% sqrt(diag(abs(ss$values))))

  ## LA

  au_la <- laplacian.spectral.embedding(g, no=no, which="la", type="I-DAD",
                                        scaled=TRUE)
  as_la <- laplacian.spectral.embedding(g, no=no, which="la", type="I-DAD",
                                        scaled=FALSE)

  expect_that(au_la$D, equals(abs(D[1:no])))
  expect_that(std(au_la$X), equals(std(X[,1:no])))
  expect_that(as_la$D, equals(D[1:no]))
  expect_that(std(as_la$X), equals(std(U[,1:no])))

  ## LM

  au_lm <- laplacian.spectral.embedding(g, no=no, which="lm", type="I-DAD",
                                        scaled=TRUE)
  as_lm <- laplacian.spectral.embedding(g, no=no, which="lm", type="I-DAD",
                                        scaled=FALSE)

  expect_that(au_lm$D, equals(mag_sort(D)[1:no]))
  expect_that(std(au_lm$X), equals(std(X[,mag_order(D)][,1:no])))
  expect_that(as_lm$D, equals(mag_sort(D)[1:no]))
  expect_that(std(as_lm$X), equals(std(U[,mag_order(D)][,1:no])))

  ## SA

  au_sa <- laplacian.spectral.embedding(g, no=no, which="sa", type="I-DAD",
                                        scaled=TRUE)
  as_sa <- laplacian.spectral.embedding(g, no=no, which="sa", type="I-DAD",
                                        scaled=FALSE)

  expect_that(au_sa$D, equals(D[vcount(g)-1:no+1]))
  expect_that(std(au_sa$X), equals(std(X[,vcount(g)-1:no+1])))
  expect_that(as_sa$D, equals(D[vcount(g)-1:no+1]))
  expect_that(std(as_sa$X), equals(std(U[,vcount(g)-1:no+1])))
})

test_that("Directed, unweighted case works", {
  library(igraph)
  g <- graph.ring(10, directed=TRUE)

  expect_that(laplacian.spectral.embedding(g, no=3),
              throws_error("not implemented for directed"))
})

test_that("Directed, weighted case works", {
  library(igraph)
  g <- graph.ring(10, directed=TRUE)
  E(g)$weight <- seq_len(ecount(g))

  expect_that(laplacian.spectral.embedding(g, no=3),
              throws_error("not implemented for directed"))

})
