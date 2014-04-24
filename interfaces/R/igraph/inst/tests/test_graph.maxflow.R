
context("graph.maxflow")

test_that("graph.maxflow works", {
  library(igraph)
  E <- rbind( c(1,3,3), c(3,4,1), c(4,2,2), c(1,5,1), c(5,6,2), c(6,2,10))
  colnames(E) <- c("from", "to", "capacity")
  g1 <- graph.data.frame(as.data.frame(E))
  fl <- graph.maxflow(g1, source="1", target="2")
  expect_that(fl$value, equals(2))
  expect_that(fl$flow, equals(rep(1, 6)))
  expect_that(sort(fl$cut), equals(c(2,4)))
  expect_that(sort(fl$partition1), equals(1:2))
  expect_that(sort(fl$partition2), equals(3:6))
})
