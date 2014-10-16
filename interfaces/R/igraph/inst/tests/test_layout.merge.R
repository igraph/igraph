
context("merge_coords")

test_that("merge_coords works", {

  library(igraph)
  set.seed(42)

  g <- list(make_ring(10), make_ring(5))
  l <- lapply(g, layout_with_mds)
  l

  lm <- merge_coords(g, l)
  expect_that(is.matrix(lm), is_true())
  expect_that(ncol(lm), equals(2))
  expect_that(nrow(lm), equals(sum(sapply(g, vcount))))

##########

  ## Stress test
  for (i in 1:10) {
    g <- sample_gnp(100, 2/100)
    l <- layout_with_mds(g)
    expect_that(dim(l), equals(c(vcount(g), 2)))
  }

})
