
context("Make API")

test_that("make_ works, order of arguments does not matter", {

  g0 <- make_undirected_graph(1:10)
  g1 <- make_(undirected_graph(1:10))
  g2 <- make_(undirected_graph(), 1:10)
  g3 <- make_(1:10, undirected_graph())

  expect_true(identical_graphs(g0, g1))
  expect_true(identical_graphs(g0, g2))
  expect_true(identical_graphs(g0, g3))
  
})

test_that("sample_, graph_ also work", {

  g0 <- make_undirected_graph(1:10)
  g1 <- sample_(undirected_graph(1:10))
  g2 <- sample_(undirected_graph(), 1:10)
  g3 <- sample_(1:10, undirected_graph())
  
  expect_true(identical_graphs(g0, g1))
  expect_true(identical_graphs(g0, g2))
  expect_true(identical_graphs(g0, g3))

  g4 <- graph_(undirected_graph(1:10))
  g5 <- graph_(undirected_graph(), 1:10)
  g6 <- graph_(1:10, undirected_graph())

  expect_true(identical_graphs(g0, g4))
  expect_true(identical_graphs(g0, g5))
  expect_true(identical_graphs(g0, g6))

})

test_that("error messages are proper", {

  expect_error(make_(), "Don't know how to make_")
  expect_error(make_(1:10), "Don't know how to make_")

  expect_error(graph_(), "Don't know how to graph_")
  expect_error(graph_(1:10), "Don't know how to graph_")
  expect_error(graph_(directed_graph(), directed_graph()),
               "Don't know how to graph_")

  expect_error(sample_(), "Don't know how to sample_")
  expect_error(sample_(1:10), "Don't know how to sample_")
  expect_error(sample_(directed_graph(), directed_graph()),
               "Don't know how to sample_")

})

test_that("we pass arguments unevaluated", {

  g0 <- graph_from_literal(A -+ B:C)
  g1 <- graph_(from_literal(A -+ B:C))
  expect_true(identical_graphs(g0, g1))

})
