
context("canonical_permutation")

test_that("canonical_permutation works", {
  library(igraph)

  g1 <- sample_gnm(10, 20)
  cp1 <- canonical_permutation(g1)
  cf1 <- permute(g1, cp1$labeling)
     
  ## Do the same with a random permutation of it
  g2 <- permute(g1, sample(vcount(g1)))
  cp2 <- canonical_permutation(g2)
  cf2 <- permute(g2, cp2$labeling)
     
  ## Check that they are the same
  el1 <- as_edgelist(cf1)
  el2 <- as_edgelist(cf2)
  el1 <- el1[ order(el1[,1], el1[,2]), ]
  el2 <- el2[ order(el2[,1], el2[,2]), ]

  expect_that(el1, equals(el2))
})
