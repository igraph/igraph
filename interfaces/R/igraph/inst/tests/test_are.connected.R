
context("is_adjacent_to")

test_that("is_adjacent_to works", {
  library(igraph)

  g <- graph_from_literal( A-B-C, B-D )
  expect_that(is_adjacent_to(g, "A", "B"), is_true())
  expect_that(is_adjacent_to(g, "B", "A"), is_true())
  expect_that(is_adjacent_to(g, "A", "D"), is_false())

  g2 <- graph( c(1,2, 2,3, 3,4), dir=FALSE )
  expect_that(is_adjacent_to(g2, 1,2), is_true())
  expect_that(is_adjacent_to(g2, 3,2), is_true())
  expect_that(is_adjacent_to(g2, 4,1), is_false())
  
  g3 <- graph_from_literal( A-+B-+C, B-+D )
  expect_that(is_adjacent_to(g3, "A", "C"), is_false())
  expect_that(is_adjacent_to(g3, "A", "B"), is_true())
  expect_that(is_adjacent_to(g3, "B", "A"), is_false())
})



