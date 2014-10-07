
context("is_bipartite")

test_that("is_bipartite works", {

  library(igraph)

  I <- matrix(sample(0:1, 35, replace=TRUE, prob=c(3,1)), nc=5)
  g <- graph_from_incidence_matrix(I)
  expect_that(bipartite_mapping(g)$res, is_true())

  set.seed(42)
  I <- matrix(sample(0:1, 35, replace=TRUE, prob=c(3,1)), nc=5)
  g <- graph_from_incidence_matrix(I)
  expect_that(bipartite_mapping(g),
              equals(list(res=TRUE, type=c(rep(FALSE, 7), rep(TRUE, 5)))))
})
