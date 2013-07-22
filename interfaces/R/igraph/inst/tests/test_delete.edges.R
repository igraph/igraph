
context("delete.edges")

test_that("delete.edges works", {
  library(igraph)
  g <- graph.formula(A:B:C - D:E:F, D-E-F)
  g2 <- delete.edges(g, E(g, P=c("D", "E")))
  expect_that(as.matrix(g2[]),
              is_equivalent_to(cbind(c(0,0,0,1,1,1), c(0,0,0,1,1,1),
                                     c(0,0,0,1,1,1), c(1,1,1,0,0,0),
                                     c(1,1,1,0,0,1), c(1,1,1,0,1,0))))
})
