
context("add_edges")

test_that("add_edges keeps edge id order", {
  library(igraph)
  g <- make_empty_graph(10)
  edges <- c(1,2, 2,3, 3,4, 1,6, 1,7, 9,10)
  g2 <- add_edges(g, edges)
  ec <- ecount(g2)
  ec2 <- length(edges)/2
  expect_equal(ec, ec2)
  expect_equal(get.edge.ids(g2, edges), seq_len(length(edges)/2))
})

test_that("add_edges adds attributes", {
  library(igraph)
  g <- make_empty_graph(10)
  g3 <- add_edges(g, (edges <- c(1,5, 2,6, 3,10, 4,5)),
                  attr=list(weight=(weights <- c(1,2,1,-1))) )
  expect_that(ecount(g3), equals(length(edges)/2))
  expect_that(get.edge.ids(g3, edges), equals(seq_len(length(edges)/2)))
  expect_that(E(g3)$weight, equals(weights))
})

test_that("add_edges unknown attributes to NA", {
  library(igraph)
  g <- make_empty_graph(10)
  g2 <- add_edges(g, (edges <- c(1,2, 2,3, 3,4, 1,6, 1,7, 9,10)) )
  g4 <- add_edges(g2, c(1,4, 4,6, 7,1), attr=list(weight=c(-1,1,-2.5)))
  expect_that(all(is.na(E(g4)$weight[seq_len(length(edges)/2)])), is_true())
})

test_that("add_edges appends attributes properly", {
  library(igraph)
  g <- make_empty_graph(10)
  g3 <- add_edges(g, (edges1 <- c(1,5, 2,6, 3,10, 4,5)),
                  attr=list(weight=(weights1 <- c(1,2,1,-1))) )
  g5 <- add_edges(g3, (edges2 <- c(10,9, 10,10, 1,1)),
                  attr=list(weight=(weights2 <- c(100,100,100))) )
  expect_that(E(g5)$weight, equals(c(weights1, weights2)))
})

test_that("add_edges signals error for zero vertex ids", {
  library(igraph)
  g <- make_full_graph(5) %du% make_full_graph(5) %du% make_full_graph(5)
  expect_that(add_edges(g, c(0,5, 0,10, 5,10)),
              throws_error("Invalid vertex id"))
})
