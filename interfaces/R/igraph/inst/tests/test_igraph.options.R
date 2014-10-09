
context("igraph_options")

test_that("igraph_options works", {

  library(igraph)

  igraph_options(verbose=TRUE)
  expect_that(igraph_opt("verbose"), is_true())
  
  igraph_options(verbose=FALSE)
  expect_that(igraph_opt("verbose"), is_false())

})
