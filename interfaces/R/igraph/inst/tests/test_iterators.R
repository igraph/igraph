
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

test_that("complex attributes work", {
  library(igraph)

  g <- graph.ring(10)
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
