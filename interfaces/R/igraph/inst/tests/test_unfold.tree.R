
context("unfold.tree")

test_that("unfold.tree works", {
  
  library(igraph)
  
  g <- graph.tree(7, 2)
  g <- add.edges(g, c(2,7, 1,4))
  g2 <- unfold.tree(g, roots=1)
  expect_that(graph.isomorphic(g2$tree, graph(c(1,2, 1,3, 2,8, 2,5, 3,6,
                                                3,9, 2,7, 1,4))), is_true())
  expect_that(g2$vertex_index, equals(c(1,2,3,4,5,6,7,4,7)))
})
