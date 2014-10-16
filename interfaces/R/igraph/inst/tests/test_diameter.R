
context("diameter")

test_that("diameter works", {

  library(igraph)

  gc <- function(graph) {
    clu <- components(graph)
    induced_subgraph(graph, which(clu$membership==which.max(clu$csize)))
  }

#### Undirected
  
  g <- gc(sample_gnp(30, 3/30))
  sp <- distances(g)
  expect_that(max(sp), equals(diameter(g)))

  g <- gc(sample_gnp(100, 1/100))
  sp <- distances(g)
  sp[sp==Inf] <- NA
  expect_that(max(sp, na.rm=TRUE), equals(diameter(g)))

#### Directed

  g <- sample_gnp(30, 3/30, dir=TRUE)
  sp <- distances(g, mode="out")
  sp[sp==Inf] <- NA
  expect_that(max(sp, na.rm=TRUE), equals(diameter(g, unconnected=TRUE)))

#### Weighted

  E(g)$weight <- sample(1:10, ecount(g), replace=TRUE)
  sp <- distances(g, mode="out")
  sp[sp==Inf] <- NA
  expect_that(max(sp, na.rm=TRUE), equals(diameter(g, unconnected=TRUE)))

#### Bug #680538

  g <- make_tree(30, mode="undirected")
  E(g)$weight <- 2
  expect_that(diameter(g, unconnected=FALSE), equals(16))
})
