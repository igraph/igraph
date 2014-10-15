
context("New layout API")

test_that("two step layouting works", {

  g <- ring(10)
  l1 <- layout_as_star(g)
  l2 <- lay_out(g, as_star())
  expect_identical(l1, l2)

})

test_that("parameters go through", {

  g <- ring(10)
  l1 <- layout_as_star(g, center = 5)
  l2 <- lay_out(g, as_star(center = 5))
  expect_identical(l1, l2)

})

test_that("parameters are evaluated early", {

  g <- ring(10)
  l1 <- layout_as_star(g, center = 5)

  cc <- 5
  spec <- as_star(center = cc)
  cc <- 10
  l2 <- lay_out(g, spec)
  expect_identical(l1, l2)

})

test_that("piping form is OK, too", {

  g <- ring(10)
  l1 <- layout_as_star(g, center = 5)
  l2 <- g %>%
    lay_out(as_star(center = 5))
  expect_identical(l1, l2)

})

test_that("add_layout works", {

  g <- ring(10)
  l1 <- layout_as_star(g, center = 5)
  l2 <- add_layout(g, as_star(center = 5))$layout
  expect_identical(l1, l2)

  l3 <- g %>%
    add_layout(as_star(center = 5)) %>%
    graph_attr("layout")
  expect_identical(l1, l3)

})
