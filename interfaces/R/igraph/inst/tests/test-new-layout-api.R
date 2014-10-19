
`%>%` <- magrittr::`%>%`

context("New layout API")

test_that("two step layouting works", {

  g <- make_ring(10)
  l1 <- layout_as_star(g)
  l2 <- layout_(g, as_star())
  expect_identical(l1, l2)

})

test_that("parameters go through", {

  g <- make_ring(10)
  l1 <- layout_as_star(g, center = 5)
  l2 <- layout_(g, as_star(center = 5))
  expect_identical(l1, l2)

})

test_that("parameters are evaluated early", {

  g <- make_ring(10)
  l1 <- layout_as_star(g, center = 5)

  cc <- 5
  spec <- as_star(center = cc)
  cc <- 10
  l2 <- layout_(g, spec)
  expect_identical(l1, l2)

})

test_that("piping form is OK, too", {

  g <- make_ring(10)
  l1 <- layout_as_star(g, center = 5)
  l2 <- g %>%
    layout_(as_star(center = 5))
  expect_identical(l1, l2)

})

test_that("add_layout_ works", {

  g <- make_ring(10)
  l1 <- layout_as_star(g, center = 5)
  l2 <- add_layout_(g, as_star(center = 5))$layout
  expect_identical(l1, l2)

  l3 <- g %>%
    add_layout_(as_star(center = 5)) %>%
    graph_attr("layout")
  expect_identical(l1, l3)

})
