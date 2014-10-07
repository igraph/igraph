
context("largest_cliques")

test_that("largest_cliques works", {

  library(igraph)

  g <- sample_gnp(50,20/50)
  lc <- largest_cliques(g)

  ## TODO: this only checks that these are cliques
  expect_that(unique(sapply(lc, function(x)
                            density(induced_subgraph(g, x)))),
              equals(1))

})
