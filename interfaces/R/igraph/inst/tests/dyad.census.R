
context("dyad.census")

test_that("dyad.census works", {

  library(igraph)

  g1 <- graph.ring(10)
  expect_that(dc1 <- dyad.census(g1), gives_warning("undirected"))
  expect_that(dc1, equals(list(mut=10, asym=0, null=35)))

  g2 <- graph.ring(10, directed=TRUE, mutual=TRUE)
  dc2 <- dyad.census(g2)
  expect_that(dc2, equals(list(mut=10, asym=0, null=35)))

  g3 <- graph.ring(10, directed=TRUE, mutual=FALSE)
  dc3 <- dyad.census(g3)
  expect_that(dc3, equals(list(mut=0, asym=10, null=35)))

  g4 <- graph.empty(2000000)
  expect_that(dc4 <- dyad.census(g4), gives_warning("Integer overflow"))
  expect_that(dc4, equals(list(mut=0, asym=0, null=0)))
  
})
