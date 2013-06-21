
context("igraph.options")

test_that("igraph.options works", {

  library(igraph)

  igraph.options(verbose=TRUE)
  expect_that(getIgraphOpt("verbose"), is_true())
  
  igraph.options(verbose=FALSE)
  expect_that(getIgraphOpt("verbose"), is_false())

})
