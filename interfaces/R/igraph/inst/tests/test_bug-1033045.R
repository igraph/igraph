
context("Bug 1033045")

test_that("Minimal s-t separators work", {
  library(igraph)
  g <- graph_from_literal(a -- 1:3 -- 5 -- 2:4 -- b, 1 -- 2, 3 -- 4)
  stsep <- min_st_separators(g)
  ims <- sapply(stsep, is_min_separator, graph=g)
  expect_that(ims, equals(rep(TRUE, 9)))
})
