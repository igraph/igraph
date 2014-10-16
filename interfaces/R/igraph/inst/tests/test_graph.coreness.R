
context("coreness")

test_that("coreness works", {
  library(igraph)
  g <- make_ring(10)
  g <- add_edges(g, c(1,2, 2,3, 1,3))
  gc <- coreness(g)               
  expect_that(gc, equals(c(3,3,3,2,2,2,2,2,2,2)))
})
