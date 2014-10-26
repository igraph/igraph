
context("Data version and conversion")

test_that("we create graphs of the current version", {

  g <- make_ring(10)
  v1 <- graph_version(g)
  v2 <- graph_version()
  expect_equal(v1, v2)

})

test_that("we can upgrade from 0.4.0 to 0.8.0", {

  g <- make_ring(10)
  g <- unclass(g)
  g[[10]] <- NULL
  class(g) <- "igraph"

  expect_equal(graph_version(g), "0.4.0")

  g2 <- upgrade_graph(g)
  expect_equal(graph_version(g2), "0.8.0")

})
