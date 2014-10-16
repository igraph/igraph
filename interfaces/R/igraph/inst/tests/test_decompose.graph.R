
context("decompose")

test_that("decompose works", {
  library(igraph)
  g <- sample_gnp(1000, 1/1500)
  G <- decompose(g)
  clu <- components(g)
  Gsizes <- sapply(G, vcount)
  expect_that(sort(clu$csize), equals(sort(Gsizes)))
})

test_that("decompose works for many components", {
  library(igraph)
  g <- make_empty_graph(50001)
  tmp <- decompose(g)
  expect_that(1, equals(1))
})

test_that("decompose works for many components and attributes", {
  library(igraph)
  g <- make_empty_graph(50001)
  V(g)$name <- 1:vcount(g)
  tmp <- decompose(g)
  expect_that(1, equals(1))
})

test_that("decompose keeps attributes", {
  library(igraph)
  g <- make_ring(10) + make_ring(5)
  V(g)$name <- letters[1:(10+5)]
  E(g)$name <- apply(as_edgelist(g), 1, paste, collapse="-")
  d <- decompose(g)
  d <- d[order(sapply(d, vcount))]

  expect_that(length(d), equals(2))
  expect_that(sapply(d, vcount), equals(c(5,10)))
  expect_that(V(d[[1]])$name, equals(letters[1:5+10]))
  expect_that(V(d[[2]])$name, equals(letters[1:10]))
  e1 <- apply(as_edgelist(d[[1]]), 1, paste, collapse="-")
  e2 <- apply(as_edgelist(d[[2]]), 1, paste, collapse="-")
  expect_that(E(d[[1]])$name, equals(e1))
  expect_that(E(d[[2]])$name, equals(e2))
})
