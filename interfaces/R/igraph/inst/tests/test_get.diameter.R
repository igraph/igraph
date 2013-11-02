
context("get.diameter")

test_that("get.diameter works", {

  library(igraph)

  g <- graph.ring(10)
  E(g)$weight <- sample(seq_len(ecount(g)))
  d <- diameter(g)
  gd <- get.diameter(g)
  sp <- shortest.paths(g)

  expect_that(d, equals(max(sp)))
  expect_that(sp[ gd[1], gd[length(gd)] ], equals(d))

  d <- diameter(g, weights=NA)
  gd <- get.diameter(g, weights=NA)
  sp <- shortest.paths(g, weights=NA)
  
  expect_that(d, equals(max(sp)))
  length(gd) == d + 1
  expect_that(sp[ gd[1], gd[length(gd)] ], equals(d))
})

