
context("make_graph")

test_that("make_graph works", {

  g <- make_graph(1:10)
  g2 <- make_empty_graph(n = 10) + edges(1:10)
  expect_true(identical_graphs(g, g2))

})

test_that("make_graph works for numeric edges and isolates", {

  g <- make_graph(1:10, n = 20)
  g2 <- make_empty_graph(n = 20) + edges(1:10)
  expect_true(identical_graphs(g, g2))

})

test_that("make_graph handles names", {

  g <- make_graph(letters[1:10])
  g2 <- make_empty_graph() + vertices(letters[1:10]) + edges(letters[1:10])
  expect_true(identical_graphs(g, g2))

})

test_that("make_graph hadles names and isolates", {

  g <- make_graph(letters[1:10], isolates = letters[11:20])
  g2 <- make_empty_graph() + vertices(letters[1:20]) + edges(letters[1:10])
  expect_true(identical_graphs(g, g2))

})

test_that("make_graph gives warning for ignored arguments", {

  expect_warning(
    make_graph(letters[1:10], n = 10)
  )

  expect_warning(
    make_graph(1:10, isolates = 11:12)
  )

})
