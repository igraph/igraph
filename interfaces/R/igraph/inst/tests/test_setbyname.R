
context("operators by name")

## TODO: we only test thay they run, but not correctness!

test_that("set operations work by name", {

  library(igraph)

  ## undirected

  g1 <- graph.ring(10)
  V(g1)$name <- rev(letters[1:10])

  g2 <- graph.star(11, center=11, mode="undirected") %u% graph.ring(10)
  V(g2)$name <- letters[1:11]

  graph.union.by.name(g1, g2)
  graph.intersection.by.name(g1, g2)
  graph.difference.by.name(g2, g1)

  ## directed

  g1 <- graph.ring(10, directed=TRUE)
  V(g1)$name <- rev(letters[1:10])

  g2 <- graph.star(11, center=11) %u% graph.ring(10, directed=TRUE,
                         mutual=TRUE)
  V(g2)$name <- letters[1:11]

  graph.union.by.name(g1, g2)
  graph.intersection.by.name(g1, g2)
  graph.difference.by.name(g2, g1)

  ## empty graph

  g1 <- graph.ring(10, directed=TRUE)
  V(g1)$name <- letters[1:10]

  g2 <- graph.empty(5, directed=TRUE)
  V(g2)$name <- letters[1:5]

  graph.union.by.name(g1, g2)
  graph.intersection.by.name(g1, g2)
  graph.difference.by.name(g1, g2)

})
