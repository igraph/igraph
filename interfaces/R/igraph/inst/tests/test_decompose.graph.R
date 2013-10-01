
context("decompose.graph")

test_that("decompose.graph works", {
  library(igraph)
  g <- erdos.renyi.game(1000, 1/1500)
  G <- decompose.graph(g)
  clu <- clusters(g)
  Gsizes <- sapply(G, vcount)
  expect_that(sort(clu$csize), equals(sort(Gsizes)))
})

test_that("decompose.graph works for many components", {
  library(igraph)
  g <- graph.empty(50001)
  tmp <- decompose.graph(g)
  expect_that(1, equals(1))
})

test_that("decompose.graph works for many components and attributes", {
  library(igraph)
  g <- graph.empty(50001)
  V(g)$name <- 1:vcount(g)
  tmp <- decompose.graph(g)
  expect_that(1, equals(1))
})

test_that("decompose.graph keeps attributes", {
  library(igraph)
  g <- graph.ring(10) + graph.ring(5)
  V(g)$name <- letters[1:(10+5)]
  E(g)$name <- apply(get.edgelist(g), 1, paste, collapse="-")
  d <- decompose.graph(g)
  d <- d[order(sapply(d, vcount))]

  expect_that(length(d), equals(2))
  expect_that(sapply(d, vcount), equals(c(5,10)))
  expect_that(V(d[[1]])$name, equals(letters[1:5+10]))
  expect_that(V(d[[2]])$name, equals(letters[1:10]))
  e1 <- apply(get.edgelist(d[[1]]), 1, paste, collapse="-")
  e2 <- apply(get.edgelist(d[[2]]), 1, paste, collapse="-")
  expect_that(E(d[[1]])$name, equals(e1))
  expect_that(E(d[[2]])$name, equals(e2))
})
