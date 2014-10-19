
context("cluster_louvain")

test_that("cluster_louvain works", {

  library(igraph)

  g <- make_graph("Zachary")
  mc <- cluster_louvain(g)
  
  expect_that(as.vector(membership(mc)),
              equals(c(2, 2, 2, 2, 1, 1, 1, 2, 4, 2, 1, 2, 2, 2, 4, 4,
                       1, 2, 4, 2, 4, 2, 4, 3, 3, 3, 4, 3, 3, 4, 4, 3,
                       4, 4) ))
  expect_that(modularity(g, mc$membership), equals(max(mc$modularity)))
  expect_that(length(mc), equals(4))
  expect_that(sizes(mc),
              equals(structure(c(5L, 12L, 6L, 11L), .Dim = 4L,
                .Dimnames = structure(list(`Community sizes` = c("1",
                     "2", "3", "4")), .Names = "Community sizes"),
                               class = "table") ))
})
