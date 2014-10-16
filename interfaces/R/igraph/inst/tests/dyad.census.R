
context("dyad_census")

test_that("dyad_census works", {

  library(igraph)

  g1 <- make_ring(10)
  expect_that(dc1 <- dyad_census(g1), gives_warning("undirected"))
  expect_that(dc1, equals(list(mut=10, asym=0, null=35)))

  g2 <- make_ring(10, directed=TRUE, mutual=TRUE)
  dc2 <- dyad_census(g2)
  expect_that(dc2, equals(list(mut=10, asym=0, null=35)))

  g3 <- make_ring(10, directed=TRUE, mutual=FALSE)
  dc3 <- dyad_census(g3)
  expect_that(dc3, equals(list(mut=0, asym=10, null=35)))

  g4 <- make_empty_graph(2000000)
  expect_that(dc4 <- dyad_census(g4), gives_warning("Integer overflow"))
  expect_that(dc4, equals(list(mut=0, asym=0, null=0)))
  
})
