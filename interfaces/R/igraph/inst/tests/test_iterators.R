
context("iterators")

test_that("iterators work", {

  library(igraph)

  ## Create a small ring graph, assign attributes
  ring <- graph.formula( A-B-C-D-E-F-G-A )
  E(ring)$weight <- seq_len(ecount(ring))

  ## Selection based on attributes
  expect_that(sort(E(ring)[ weight < 4 ]$weight), equals(1:3))
  expect_that(V(ring)[ c("A", "C") ]$name, equals(c("A", "C")))

  ## TODO: %--%, %->%, other special functions

})
