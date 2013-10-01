
context("graph.isoclass")

test_that("graph.isoclass works", {

  library(igraph)

  g1 <- graph.isocreate(3, 10)
  g2 <- graph.isocreate(3, 11)
  expect_that(graph.isoclass(g1), equals(10))
  expect_that(graph.isoclass(g2), equals(11))

  g1 <- add.vertices(g1, 3)
  expect_that(graph.isoclass.subgraph(g1, 1:3), equals(10))
  expect_that(graph.isoclass.subgraph(g1 %du% g2, 1:3), equals(10))
  expect_that(graph.isoclass.subgraph(g1 %du% g2, 7:9), equals(11))

})
