
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

test_that("brackering works", {
  library(igraph)

  g <- graph(c(1,2, 1,3, 3,4))
  g <- set.vertex.attribute(g, name="weight", value=1:vcount(g))
  g <- set.edge.attribute(g, name="weight", value=1:ecount(g))
  g <- set.graph.attribute(g, name="name", "foo")

  graph2 <- set.vertex.attribute(g, name="weight",
                                 value=rep(1, vcount(g)))
  graph2 <- set.edge.attribute(g, name="weight",
                               value=rep(1, ecount(g)))
  graph2 <- set.graph.attribute(g, name="name", "foobar")

  expect_that(get.vertex.attribute(g, name="weight"),
              equals(1:4))
  expect_that(get.edge.attribute(g, name="weight"),
              equals(1:3))
  expect_that(get.graph.attribute(g, name="name"), equals("foo"))
})

test_that("brackering works with a function", {
  library(igraph)
  library(testthat)

  g <- graph(c(1,2, 1,3, 3,4))
  g <- set.vertex.attribute(g, name="weight", value=1:vcount(g))
  g <- set.edge.attribute(g, name="weight", value=1:ecount(g))
  g <- set.graph.attribute(g, name="name", "foo")

  run.test <- function(graph) {
    graph2 <- set.vertex.attribute(graph, name="weight",
                                   value=rep(1, vcount(graph)))
    graph2 <- set.edge.attribute(graph, name="weight",
                                   value=rep(1, ecount(graph)))
    graph2 <- set.graph.attribute(graph, name="name", "foobar")
  }

  g2 <- run.test(g)
  expect_that(get.vertex.attribute(g, name="weight"),
              equals(1:4))
  expect_that(get.edge.attribute(g, name="weight"),
              equals(1:3))
  expect_that(get.graph.attribute(g, name="name"), equals("foo"))
})

test_that("brackering works with shortcuts", {
  library(igraph)

  g <- graph(c(1,2, 1,3, 3,4))
  g <- set.vertex.attribute(g, name="weight", value=1:vcount(g))
  g <- set.edge.attribute(g, name="weight", value=1:ecount(g))
  g <- set.graph.attribute(g, name="name", "foo")

  run.test <- function(graph) {
    V(graph)$weight <- rep(1, vcount(graph))
    E(graph)$weight <- rep(1, ecount(graph))
    graph$name <- "foobar"
  }

  g2 <- run.test(g)
  expect_that(get.vertex.attribute(g, name="weight"),
              equals(1:4))
  expect_that(get.edge.attribute(g, name="weight"),
              equals(1:3))
  expect_that(get.graph.attribute(g, name="name"), equals("foo"))
})

## TODO: subsetting


