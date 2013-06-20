
context("layout.mds")

test_that("layout.mds works", {

  library(igraph)

  ## A tree

  g <- graph.tree(10, 2, "undirected")

  mymds <- function(g) { 
    sp <- shortest.paths(g)
    sp <- sp * sp
    sp <- sp - rowMeans(sp) - rep(rowMeans(sp), each=nrow(sp)) + mean(sp)
    sp <- sp / -2
    ei <- eigen(sp)
    va <- sqrt(abs(ei$values[1:2]))
    ei$vectors[,1:2] * rep(va, each=nrow(sp))
  }

  expect_that(mymds(g), equals(layout.mds(g)))

  ## plot(g, layout=ll)

  ## A graph with multiple components, just test that it runs

  set.seed(42)
  g <- graph.ring(10) + graph.ring(3)
  expect_that(ncol(layout.mds(g)), equals(2))
  
  ## Small stress test

  for (i in 1:10) {
    g <- erdos.renyi.game(100, 2/100)
    l <- layout.mds(g)
    expect_that(ncol(l), equals(2))
  }

})
