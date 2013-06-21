
context("canonical.permutation")

test_that("canonical.permutation works", {
  library(igraph)

  g1 <- erdos.renyi.game(10, 20, type="gnm")
  cp1 <- canonical.permutation(g1)
  cf1 <- permute.vertices(g1, cp1$labeling)
     
  ## Do the same with a random permutation of it
  g2 <- permute.vertices(g1, sample(vcount(g1)))
  cp2 <- canonical.permutation(g2)
  cf2 <- permute.vertices(g2, cp2$labeling)
     
  ## Check that they are the same
  el1 <- get.edgelist(cf1)
  el2 <- get.edgelist(cf2)
  el1 <- el1[ order(el1[,1], el1[,2]), ]
  el2 <- el2[ order(el2[,1], el2[,2]), ]

  expect_that(el1, equals(el2))
})
