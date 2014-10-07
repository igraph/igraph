
context("Triangles")

test_that("Listing triangles works", {

  library(igraph)
  g1 <- empty_graph(directed=TRUE)
  g2 <- empty_graph(directed=FALSE)
  expect_that(triangles(g1), equals(numeric()))
  expect_that(triangles(g2), equals(numeric()))

  g3 <- empty_graph(n=1, directed=TRUE)
  g4 <- empty_graph(n=1, directed=FALSE)
  expect_that(triangles(g3), equals(numeric()))
  expect_that(triangles(g4), equals(numeric()))

  g5 <- empty_graph(n=100, directed=TRUE)
  g6 <- empty_graph(n=100, directed=FALSE)
  expect_that(triangles(g5), equals(numeric()))
  expect_that(triangles(g6), equals(numeric()))
  
  g7 <- ring(3, directed=FALSE)
  g8 <- ring(3, directed=TRUE)
  g9 <- graph_from_formula(A-+B:C, B-+C)
  expect_that(sort(triangles(g7)), equals(1:3))
  expect_that(sort(triangles(g8)), equals(1:3))
  expect_that(sort(triangles(g9)), equals(1:3))

  g10 <- full_graph(5, directed=FALSE)
  g11 <- full_graph(5, directed=TRUE)
  r10 <- c(1L, 2L, 5L, 1L, 2L, 3L, 1L, 2L, 4L, 1L, 3L, 5L, 1L, 3L, 4L, 
           1L, 4L, 5L, 2L, 3L, 5L, 2L, 3L, 4L, 2L, 4L, 5L, 3L, 4L, 5L)
  r11 <- c(1L, 2L, 5L, 1L, 2L, 4L, 1L, 2L, 3L, 1L, 3L, 5L, 1L, 3L, 4L, 
           1L, 4L, 5L, 2L, 4L, 5L, 2L, 3L, 5L, 2L, 3L, 4L, 3L, 4L, 5L)
  expect_that(triangles(g10), equals(r10))
  expect_that(triangles(g11), equals(r11))

})
