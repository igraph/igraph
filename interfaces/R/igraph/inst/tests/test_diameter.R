
context("diameter")

test_that("diameter works", {

  library(igraph)

  gc <- function(graph) {
    clu <- clusters(graph)
    induced.subgraph(graph, which(clu$membership==which.max(clu$csize)))
  }

#### Undirected
  
  g <- gc(erdos.renyi.game(30, 3/30))
  sp <- shortest.paths(g)
  expect_that(max(sp), equals(diameter(g)))

  g <- gc(erdos.renyi.game(100, 1/100))
  sp <- shortest.paths(g)
  sp[sp==Inf] <- NA
  expect_that(max(sp, na.rm=TRUE), equals(diameter(g)))

#### Directed

  g <- erdos.renyi.game(30, 3/30, dir=TRUE)
  sp <- shortest.paths(g, mode="out")
  sp[sp==Inf] <- NA
  expect_that(max(sp, na.rm=TRUE), equals(diameter(g, unconnected=TRUE)))

#### Weighted

  E(g)$weight <- sample(1:10, ecount(g), replace=TRUE)
  sp <- shortest.paths(g, mode="out")
  sp[sp==Inf] <- NA
  expect_that(max(sp, na.rm=TRUE), equals(diameter(g, unconnected=TRUE)))

#### Bug #680538

  g <- graph.tree(30, mode="undirected")
  E(g)$weight <- 2
  expect_that(diameter(g, unconnected=FALSE), equals(16))
})
