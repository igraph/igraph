
context("authority.score")

test_that("authority score works", {
  library(igraph)

  ashs <- function(graph, as=TRUE) {
    mscale <- function(x) {
      if (sd(x)!=0) { x <- scale(x) }
      if (x[1] < 0) { x <- -x       }
      x
    }
    A <- get.adjacency(graph, sparse=FALSE)
    if (as) { 
      s1 <- eigen(t(A) %*% A)$vectors[,1]
      s2 <- authority.score(graph)$vector
    } else {
      s1 <- eigen(A %*% t(A))$vectors[,1]
      s2 <- hub.score(graph)$vector
    }
    expect_that(mscale(s1), is_equivalent_to(mscale(s2)))
  }

  g1 <- ba.game(100, m=10)
  ashs(g1)
  ashs(g1, as=FALSE)

  g2 <- erdos.renyi.game(100, 2/100)
  ashs(g2)
  ashs(g2, as=FALSE)
})

test_that("authority scores of a ring are all one", {
  library(igraph)
  g3 <- graph.ring(100)
  expect_that(authority.score(g3)$vector, equals(rep(1, vcount(g3))))
  expect_that(hub.score(g3)$vector, equals(rep(1, vcount(g3))))
})

test_that("authority.score survives stress test", {
  library(igraph)
  library(Matrix)
  set.seed(42)

  is.principal <- function(M, lambda) {
    expect_that(eigen(M)$values[1], equals(lambda))
  }

  is.ev <- function(M, v, lambda) {
    expect_that(as.vector(M %*% v), equals(lambda * v))
  }

  is.good <- function(M, v, lambda) {
    is.principal(M, lambda)
    is.ev(M, v, lambda)
  }

  for (i in 1:100) {
    G <- erdos.renyi.game(10, sample(1:20, 1), type="gnm")
    as <- authority.score(G)
    M <- get.adjacency(G)
    is.good(t(M) %*% M, as$vector, as$value)
  }

  for (i in 1:100) {
    G <- erdos.renyi.game(10, sample(1:20, 1), type="gnm")
    hs <- hub.score(G)
    M <- get.adjacency(G)
    is.good(M %*% t(M), hs$vector, hs$value)
  }
})
