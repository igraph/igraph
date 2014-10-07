
context("largest_ivs")

test_that("largest_ivs works", {

  library(igraph)

  g <- sample_gnp(50, 0.8)
  livs <- largest_ivs(g)
  expect_that(unique(sapply(livs, length)),
              equals(ivs_size(g)))

  ec <- sapply(seq_along(livs), function(x)
               ecount(induced_subgraph(g, livs[[x]])))
  expect_that(unique(ec), equals(0))

  ## TODO: check that they are largest
})
