
context("attributes")

test_that("assigning and querying attributes work", {
  library(igraph)

  ## Create a small ring graph, assign attributes
  ring <- graph_from_literal( A-B-C-D-E-F-G-A )
  E(ring)$weight <- seq_len(ecount(ring))

  ## Query attributes
  expect_that(V(ring)$name, equals(LETTERS[seq_len(vcount(ring))]))
  expect_that(E(ring)$weight, equals(seq_len(ecount(ring))))
})

test_that("brackering works", {
  library(igraph)

  g <- graph(c(1,2, 1,3, 3,4))
  g <- set_vertex_attr(g, name="weight", value=1:vcount(g))
  g <- set_edge_attr(g, name="weight", value=1:ecount(g))
  g <- set_graph_attr(g, name="name", "foo")

  graph2 <- set_vertex_attr(g, name="weight",
                                 value=rep(1, vcount(g)))
  graph2 <- set_edge_attr(g, name="weight",
                               value=rep(1, ecount(g)))
  graph2 <- set_graph_attr(g, name="name", "foobar")

  expect_that(vertex_attr(g, name="weight"),
              equals(1:4))
  expect_that(edge_attr(g, name="weight"),
              equals(1:3))
  expect_that(graph_attr(g, name="name"), equals("foo"))
})

test_that("brackering works with a function", {
  library(igraph)
  library(testthat)

  g <- graph(c(1,2, 1,3, 3,4))
  g <- set_vertex_attr(g, name="weight", value=1:vcount(g))
  g <- set_edge_attr(g, name="weight", value=1:ecount(g))
  g <- set_graph_attr(g, name="name", "foo")

  run.test <- function(graph) {
    graph2 <- set_vertex_attr(graph, name="weight",
                                   value=rep(1, vcount(graph)))
    graph2 <- set_edge_attr(graph, name="weight",
                                   value=rep(1, ecount(graph)))
    graph2 <- set_graph_attr(graph, name="name", "foobar")
  }

  g2 <- run.test(g)
  expect_that(vertex_attr(g, name="weight"),
              equals(1:4))
  expect_that(edge_attr(g, name="weight"),
              equals(1:3))
  expect_that(graph_attr(g, name="name"), equals("foo"))
})

test_that("brackering works with shortcuts", {
  library(igraph)

  g <- graph(c(1,2, 1,3, 3,4))
  g <- set_vertex_attr(g, name="weight", value=1:vcount(g))
  g <- set_edge_attr(g, name="weight", value=1:ecount(g))
  g <- set_graph_attr(g, name="name", "foo")

  run.test <- function(graph) {
    V(graph)$weight <- rep(1, vcount(graph))
    E(graph)$weight <- rep(1, ecount(graph))
    graph$name <- "foobar"
  }

  g2 <- run.test(g)
  expect_that(vertex_attr(g, name="weight"),
              equals(1:4))
  expect_that(edge_attr(g, name="weight"),
              equals(1:3))
  expect_that(graph_attr(g, name="name"), equals("foo"))
})

## TODO: subsetting

test_that("we can query all attributes at once", {

  g <- graph(c(1,2, 1,3, 2,4))

  expect_equal(graph_attr(g), structure(list(), .Names = character(0)))
  expect_equal(vertex_attr(g), list())
  expect_equal(edge_attr(g), list())

  g$name <- "toy"
  g$layout <- cbind(1:4, 1:4)
  V(g)$name <- letters[1:4]
  V(g)$color <- rainbow(4)
  E(g)$weight <- 1:3
  E(g)$label <- LETTERS[1:3]

  expect_equal(graph_attr(g), list(name = "toy", layout = cbind(1:4, 1:4)))
  expect_equal(vertex_attr(g), list(name = letters[1:4], color = rainbow(4)))
  expect_equal(edge_attr(g), list(weight = 1:3, label = LETTERS[1:3]))

})

test_that("we can query single attributes with the generic functions", {

  g <- graph(c(1,2, 1,3, 2,4))

  g$name <- "toy"
  g$layout <- cbind(1:4, 1:4)
  V(g)$name <- letters[1:4]
  V(g)$color <- rainbow(4)
  E(g)$weight <- 1:3
  E(g)$label <- LETTERS[1:3]

  expect_equal(graph_attr(g, "name"), "toy")
  expect_equal(graph_attr(g, "layout"), cbind(1:4, 1:4))
  expect_equal(vertex_attr(g, "name"), letters[1:4])
  expect_equal(vertex_attr(g, "color"), rainbow(4))
  expect_equal(edge_attr(g, "weight"), 1:3)
  expect_equal(edge_attr(g, "label"), LETTERS[1:3])

})

test_that("we can query a subset of vertices", {

  g <- graph(c(1,2, 1,3, 2,4))

  V(g)$name <- letters[1:4]
  V(g)$color <- as.list(rainbow(4))
  E(g)$weight <- 1:3
  E(g)$label <- as.list(LETTERS[1:3])

  expect_equal(vertex_attr(g, "name", c(1,3)), letters[c(1,3)])
  expect_equal(vertex_attr(g, "color", c("a", "c")),
               as.list(rainbow(4))[c(1,3)])
  expect_equal(edge_attr(g, "weight", 2:3), 2:3)
  expect_equal(edge_attr(g, "label", 2:3), as.list(LETTERS[1:3])[2:3])

})

test_that("we can set all attributes at once", {

  g <- graph(c(1,2, 1,3, 2,4))

  g$name <- "toy"
  g$layout <- cbind(1:4, 1:4)
  V(g)$name <- letters[1:4]
  V(g)$color <- as.list(rainbow(4))
  E(g)$weight <- 1:3
  E(g)$label <- as.list(LETTERS[1:3])

  g2 <- graph(c(2,1, 3,1, 4,1))

  graph_attr(g2) <- graph_attr(g)
  expect_equal(graph_attr(g2), graph_attr(g))

  vertex_attr(g2) <- vertex_attr(g)
  expect_equal(vertex_attr(g2), vertex_attr(g))

  edge_attr(g2) <- edge_attr(g)
  expect_equal(edge_attr(g2), edge_attr(g))
})

test_that("we can set all attributes some vertices/edges", {

  g <- graph(c(1,2, 1,3, 2,4))

  V(g)$name <- letters[1:4]
  V(g)$color <- as.list(rainbow(4))
  E(g)$weight <- 1:3
  E(g)$label <- as.list(LETTERS[1:3])

  g2 <- graph(c(2,1, 3,1, 4,1, 2,5, 3,6))

  vertex_attr(g2, index = c(1, 2, 4, 5)) <- vertex_attr(g)
  expect_equal(vertex_attr(g2), list(name = c("a", "b", NA_character_,
    "c", "d", NA_character_),color = list(rainbow(4)[1], rainbow(4)[2], NULL,
    rainbow(4)[3], rainbow(4)[4], NULL)))

  edge_attr(g2, index = c(1, 3, 5)) <- edge_attr(g)
  expect_equal(edge_attr(g2), list(weight = c(1L, NA_integer_, 2L,
    NA_integer_, 3L), label = list("A", NULL, "B", NULL, "C")))

})

test_that("cannot use vs/es from another graph", {

  g <- graph.ring(10)
  g2 <- g + 1
  v <- V(g)[1:4]
  expect_error(g2 - v, "Cannot use a vertex sequence from another graph")

  e <- E(g)[1:2]
  expect_error(g2 - e, "Cannot use an edge sequence from another graph")
})
