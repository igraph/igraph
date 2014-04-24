
context("graph.bipartite")

test_that("graph.bipartite works", {
  library(igraph)

  I <- matrix(sample(0:1, 35, replace=TRUE, prob=c(3,1)), nc=5)
  g <- graph.incidence(I)

  edges <- unlist(sapply(seq_len(nrow(I)), function(x) {
    w <- which(I[x,] != 0) + nrow(I)
    if (length(w)!=0) { as.vector(rbind(x, w)) } else { numeric() }
  }))
  g2 <- graph.bipartite(seq_len(nrow(I)+ncol(I)) > nrow(I), edges)
  I2 <- get.incidence(g2)
  
  expect_that(I2, is_equivalent_to(I))
})
