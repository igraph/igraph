
context("Notable graphs")

test_that("notable graphs work with make_graph", {

  g <- make_graph("Levi")
  g2 <- graph.famous("Levi")
  expect_true(identical_graphs(g, g2))

})

test_that("make_graph for notable graphs is case insensitive", {

  g <- make_graph("Levi")
  g2 <- make_graph("levi")
  expect_true(identical_graphs(g, g2))

})

test_that("spaces are replaced in make_graph for notable graphs", {

  g <- make_graph("Krackhardt_Kite")
  g2 <- make_graph("Krackhardt kite")
  expect_true(identical_graphs(g, g2))

})

test_that("warnings are given for extra arguments in make_graph for notables", {

  g0 <- make_graph("Levi")
  expect_warning(g1 <- make_graph("Levi", n = 10))
  expect_warning(g2 <- make_graph("Levi", isolates = "foo"))
  expect_warning(g3 <- make_graph("Levi", directed = FALSE))
  expect_true(identical_graphs(g0, g1))
  expect_true(identical_graphs(g0, g2))
  expect_true(identical_graphs(g0, g3))

})
