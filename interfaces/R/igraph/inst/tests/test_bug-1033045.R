
context("Bug 1033045")

test_that("Minimal s-t separators work", {
  library(igraph)
  g <- graph.formula(a -- 1:3 -- 5 -- 2:4 -- b, 1 -- 2, 3 -- 4)
  stsep <- minimal.st.separators(g)
  ims <- sapply(stsep, is.minimal.separator, graph=g)
  expect_that(ims, equals(rep(TRUE, 9)))
})
