
context("graph.mincut")

test_that("graph.mincut works", {

  library(igraph)

  g2 <- graph( c(1,2,2,3,3,4, 1,6,6,5,5,4, 4,1) )
  E(g2)$capacity <- c(3,1,2, 10,1,3, 2)
  mc <- graph.mincut(g2, value.only=FALSE)

  expect_that(mc$value, equals(1))
  expect_that(mc$cut, equals(2))
  expect_that(mc$partition1, equals(2))
  expect_that(mc$partition2, equals(c(1,3:6)))

})
