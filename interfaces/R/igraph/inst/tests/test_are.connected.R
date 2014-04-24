
context("are.connected")

test_that("are.connected works", {
  library(igraph)

  g <- graph.formula( A-B-C, B-D )
  expect_that(are.connected(g, "A", "B"), is_true())
  expect_that(are.connected(g, "B", "A"), is_true())
  expect_that(are.connected(g, "A", "D"), is_false())

  g2 <- graph( c(1,2, 2,3, 3,4), dir=FALSE )
  expect_that(are.connected(g2, 1,2), is_true())
  expect_that(are.connected(g2, 3,2), is_true())
  expect_that(are.connected(g2, 4,1), is_false())
  
  g3 <- graph.formula( A-+B-+C, B-+D )
  expect_that(are.connected(g3, "A", "C"), is_false())
  expect_that(are.connected(g3, "A", "B"), is_true())
  expect_that(are.connected(g3, "B", "A"), is_false())
})



