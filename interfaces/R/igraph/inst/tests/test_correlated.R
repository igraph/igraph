
context("Correlated E-R random graphs")

## Not very meaningful tests. They good for testing that the
## functions run, but not much more

test_that("sample_correlated_gnp works", {

  library(igraph)
  set.seed(42)

  g <- erdos.renyi.game(10, .1)
  g2 <- sample_correlated_gnp(g, corr=1, p=g$p, perm=NULL)
  expect_that(g[], equals(g2[]))

  g3 <- sample_correlated_gnp(g, corr=0, p=g$p, perm=NULL)
  c3 <- cor(as.vector(g[]), as.vector(g3[]))
  expect_that(abs(c3) < .3, is_true())

})

test_that("sample_correlated_gnp_pair works", {

  library(igraph)
  set.seed(42)

  gp <- sample_correlated_gnp_pair(10, corr=.95, p=.1, perm=NULL)
  expect_that(abs(ecount(gp[[1]]) - ecount(gp[[2]])) < 3, is_true())

})

## Some corner cases

test_that("sample_correlated_gnp corner cases work", {

  library(igraph)
  set.seed(42)

  is.full <- function(g) {
    g2 <- graph.full(vcount(g), directed=is.directed(g))
    graph.isomorphic(g, g2)
  }

  g <- erdos.renyi.game(10, .3)
  g2 <- sample_correlated_gnp(g, corr=0.000001, p=.99999999)
  expect_that(is.full(g2), is_true())

  g3 <- sample_correlated_gnp(g, corr=0.000001, p=0.0000001)
  expect_that(ecount(g3), equals(0))
  expect_that(vcount(g3), equals(10))

  gg <- erdos.renyi.game(10, .3, directed=TRUE)
  gg2 <- sample_correlated_gnp(gg, corr=0.000001, p=.99999999)
  expect_that(is.full(gg2), is_true())

  gg3 <- sample_correlated_gnp(gg, corr=0.000001, p=0.0000001)
  expect_that(ecount(gg3), equals(0))
  expect_that(vcount(gg3), equals(10))

})

test_that("permutation works for sample_correlated_gnp", {

  library(igraph)
  set.seed(42)

  g <- erdos.renyi.game(10, .3)
  perm <- sample(vcount(g))
  g2 <- sample_correlated_gnp(g, corr=.99999, p=.3, permutation=perm)
  g <- permute.vertices(g, perm)
  expect_that(g[], equals(g2[]))

  g <- erdos.renyi.game(10, .3)
  perm <- sample(vcount(g))
  g2 <- sample_correlated_gnp(g, corr=1, p=.3, permutation=perm)
  g <- permute.vertices(g, perm)
  expect_that(g[], equals(g2[]))

})
