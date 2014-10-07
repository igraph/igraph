
context("iso_class")

test_that("iso_class works", {

  library(igraph)

  g1 <- graph_from_iso_class(3, 10)
  g2 <- graph_from_iso_class(3, 11)
  expect_that(iso_class(g1), equals(10))
  expect_that(iso_class(g2), equals(11))

  g1 <- add_vertices(g1, 3)
  expect_that(graph.isoclass.subgraph(g1, 1:3), equals(10))
  expect_that(graph.isoclass.subgraph(g1 %du% g2, 1:3), equals(10))
  expect_that(graph.isoclass.subgraph(g1 %du% g2, 7:9), equals(11))

})
