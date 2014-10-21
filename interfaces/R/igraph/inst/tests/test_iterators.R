
context("iterators")

test_that("iterators work", {

  library(igraph)

  ## Create a small ring graph, assign attributes
  ring <- graph_from_literal( A-B-C-D-E-F-G-A )
  E(ring)$weight <- seq_len(ecount(ring))

  ## Selection based on attributes
  expect_that(sort(E(ring)[ weight < 4 ]$weight), equals(1:3))
  expect_that(V(ring)[ c("A", "C") ]$name, equals(c("A", "C")))

  ## TODO: %--%, %->%, other special functions

})

test_that("complex attributes work", {
  library(igraph)

  g <- make_ring(10)
  foo <- lapply(1:vcount(g), seq, from=1)
  V(g)$foo <- foo

  V(g)$foo[[5]][1:3] <- 0
  expect_that(V(g)[(1:10)[-5]]$foo, equals(foo[-5]))
  expect_that(V(g)[[5]]$foo, equals(c(0,0,0,4,5)))
  expect_that(V(g)[5]$foo, equals(list(c(0,0,0,4,5))))

  V(g)$foo <- foo
  V(g)[[5]]$foo[1:3] <- 0
  expect_that(V(g)[(1:10)[-5]]$foo, equals(foo[-5]))
  expect_that(V(g)[[5]]$foo, equals(c(0,0,0,4,5)))
  expect_that(V(g)[5]$foo, equals(list(c(0,0,0,4,5))))

  V(g)$foo <- foo
  V(g)[5]$foo[[1]][1:3] <- 0
  expect_that(V(g)[(1:10)[-5]]$foo, equals(foo[-5]))
  expect_that(V(g)[[5]]$foo, equals(c(0,0,0,4,5)))
  expect_that(V(g)[5]$foo, equals(list(c(0,0,0,4,5))))

})

test_that("we got rid of confusing indexing by numbers", {

  g <- make_ring(10)
  V(g)$name <- letters[1:10]
  E(g)$weight <- seq(ecount(g))

  expect_equal(as.vector(V(g)[6:10][1:5]), 6:10)
  expect_equal(as.vector(E(g)[6:10][1:5]), 6:10)

})

test_that("selecting edges using vertex names works", {

  g <- make_ring(10)
  V(g)$name <- letters[1:10]

  e1 <- E(g)[c('a|b', 'c|d')]
  expect_equal(as.vector(e1), c(1,3))

})

test_that("indexing with characters work as expected", {

  g <- make_ring(10)
  V(g)$name <- letters[1:10]
  E(g)$name <- LETTERS[1:10]

  expect_equal(as.vector(V(g)[letters[3:6]]), 3:6)
  expect_equal(as.vector(E(g)[LETTERS[4:7]]), 4:7)
  ## expect_equal(as.vector(E(g)[c('a|b', 'c|d')]), c(1,3))

  expect_error(V(g)[1:5]['h'], 'Unknown vertex selected')
  expect_error(E(g)[1:5]['H'], 'Unknown edge selected')
  expect_error(E(g)[6:9]['a|b'], 'Unknown edge selected')
})
