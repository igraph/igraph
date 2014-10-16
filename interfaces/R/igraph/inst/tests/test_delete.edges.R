
context("delete_edges")

test_that("delete_edges works", {
  library(igraph)
  g <- graph_from_literal(A:B:C - D:E:F, D-E-F)
  g2 <- delete_edges(g, E(g, P=c("D", "E")))
  expect_that(as.matrix(g2[]),
              is_equivalent_to(cbind(c(0,0,0,1,1,1), c(0,0,0,1,1,1),
                                     c(0,0,0,1,1,1), c(1,1,1,0,0,0),
                                     c(1,1,1,0,0,1), c(1,1,1,0,1,0))))
})
