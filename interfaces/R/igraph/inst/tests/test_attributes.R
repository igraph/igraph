
context("attributes")

test_that("assigning and querying attributes work", {
  library(igraph)

  ## Create a small ring graph, assign attributes
  ring <- graph.formula( A-B-C-D-E-F-G-A )
  E(ring)$weight <- seq_len(ecount(ring))
  
  ## Query attributes
  expect_that(V(ring)$name, equals(LETTERS[seq_len(vcount(ring))]))
  expect_that(E(ring)$weight, equals(seq_len(ecount(ring))))
})

## TODO: subsetting


