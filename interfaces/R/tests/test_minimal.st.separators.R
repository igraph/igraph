
context("minimal.st.separators")

test_that("minimal.st.separators works", {

  library(igraph)
  g <- graph.famous("Zachary")
  msts <- minimal.st.separators(g)
  is <- sapply(msts, is.separator, graph=g)
  expect_that(unique(is), equals(TRUE))

  ## TODO: check that it is minimal
})
