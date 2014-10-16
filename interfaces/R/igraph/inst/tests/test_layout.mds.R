
context("layout_with_mds")

test_that("layout_with_mds works", {

  library(igraph)

  ## A tree

  g <- make_tree(10, 2, "undirected")

  mymds <- function(g) { 
    sp <- distances(g)
    sp <- sp * sp
    sp <- sp - rowMeans(sp) - rep(rowMeans(sp), each=nrow(sp)) + mean(sp)
    sp <- sp / -2
    ei <- eigen(sp)
    va <- sqrt(abs(ei$values[1:2]))
    ei$vectors[,1:2] * rep(va, each=nrow(sp))
  }

  expect_that(mymds(g), equals(layout_with_mds(g)))

  ## plot(g, layout=ll)

  ## A graph with multiple components, just test that it runs

  set.seed(42)
  g <- make_ring(10) + make_ring(3)
  expect_that(ncol(layout_with_mds(g)), equals(2))
  
  ## Small stress test

  for (i in 1:10) {
    g <- sample_gnp(100, 2/100)
    l <- layout_with_mds(g)
    expect_that(ncol(l), equals(2))
  }

})
