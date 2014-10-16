
context("get_diameter")

test_that("get_diameter works", {

  library(igraph)

  g <- make_ring(10)
  E(g)$weight <- sample(seq_len(ecount(g)))
  d <- diameter(g)
  gd <- get_diameter(g)
  sp <- distances(g)

  expect_that(d, equals(max(sp)))
  expect_that(sp[ gd[1], gd[length(gd)] ], equals(d))

  d <- diameter(g, weights=NA)
  gd <- get_diameter(g, weights=NA)
  sp <- distances(g, weights=NA)
  
  expect_that(d, equals(max(sp)))
  length(gd) == d + 1
  expect_that(sp[ gd[1], gd[length(gd)] ], equals(d))
})

