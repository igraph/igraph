
context("clusters")

test_that("clusters works", {
  library(igraph)
  set.seed(42)
  
  gc <- function(graph) {
    cl <- clusters(graph)
    induced.subgraph(graph, which(cl$membership==which.max(cl$csize)))
  }
  
  rg <- function(n) {
    gc(erdos.renyi.game(n, 1/n))
  }
  
  G <- lapply(1:30, function(x) rg(sample(100, 1)))
  Gsize <- sapply(G, vcount)

  allg <- graph.disjoint.union(G)
  clu <- clusters(allg)

  expect_that(as.numeric(table(clu$membership)), equals(clu$csize))
  expect_that(sort(clu$csize), equals(sort(Gsize)))
  expect_that(clu$no, equals(length(G)))
})

