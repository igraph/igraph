
context("Prefix sum tree")

test_that("Prefix sum tree works", {

  library(igraph)
  set.seed(42)
  mysample <- function(x, size, prob=NULL) {
    if (!is.null(prob)) { prob <- as.numeric(prob) }
    .Call("R_igraph_psumtree_draw", as.integer(x), as.integer(size),
          prob, PACKAGE="igraph")
  }
  S <- mysample(100, 10000)
  expect_that(range(table(S)), equals(c(69, 129)))

  S2 <- mysample(100, 10000, rep(1:2, each=50))
  expect_that(range(table(S2)[1:50]), equals(c(45, 85)))
  expect_that(range(table(S2)[51:100]), equals(c(103, 160)))

})
