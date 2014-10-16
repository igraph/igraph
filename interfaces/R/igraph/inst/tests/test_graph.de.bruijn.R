
context("make_de_bruijn_graph")

test_that("make_de_bruijn_graph works", {

  library(igraph)
  g <- make_de_bruijn_graph(2,1)
  g2 <- make_de_bruijn_graph(2,2)
  g3 <- make_line_graph(g)

  expect_that(graph.isomorphic(g3, graph(c(1,1,3,1,1,2,3,2,2,3,
                                           4,3,2,4,4,4))), is_true())
  expect_that(graph.isomorphic(g2, g3), is_true())
  
})
