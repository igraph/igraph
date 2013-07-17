
context("add.edges")

test_that("add.edges keeps edge id order", {
  library(igraph)
  g <- graph.empty(10)
  g2 <- add.edges(g, (edges <- c(1,2, 2,3, 3,4, 1,6, 1,7, 9,10)) )
  expect_that(ecount(g2), equals(length(edges)/2))
  expect_that(get.edge.ids(g2, edges), equals(seq_len(length(edges)/2)))
})

test_that("add.edges adds attributes", {
  library(igraph)
  g <- graph.empty(10)
  g3 <- add.edges(g, (edges <- c(1,5, 2,6, 3,10, 4,5)),
                  attr=list(weight=(weights <- c(1,2,1,-1))) )
  expect_that(ecount(g3), equals(length(edges)/2))
  expect_that(get.edge.ids(g3, edges), equals(seq_len(length(edges)/2)))
  expect_that(E(g3)$weight, equals(weights))
})

test_that("add.edges unknown attributes to NA", {
  library(igraph)
  g <- graph.empty(10)
  g2 <- add.edges(g, (edges <- c(1,2, 2,3, 3,4, 1,6, 1,7, 9,10)) )
  g4 <- add.edges(g2, c(1,4, 4,6, 7,1), attr=list(weight=c(-1,1,-2.5)))
  expect_that(all(is.na(E(g4)$weight[seq_len(length(edges)/2)])), is_true())
})

test_that("add.edges appends attributes properly", {
  library(igraph)
  g <- graph.empty(10)
  g3 <- add.edges(g, (edges1 <- c(1,5, 2,6, 3,10, 4,5)),
                  attr=list(weight=(weights1 <- c(1,2,1,-1))) )
  g5 <- add.edges(g3, (edges2 <- c(10,9, 10,10, 1,1)),
                  attr=list(weight=(weights2 <- c(100,100,100))) )
  expect_that(E(g5)$weight, equals(c(weights1, weights2)))
})

test_that("add.edges signals error for zero vertex ids", {
  library(igraph)
  g <- graph.full(5) %du% graph.full(5) %du% graph.full(5)
  expect_that(add.edges(g, c(0,5, 0,10, 5,10)),
              throws_error("Invalid vertex id"))
})
