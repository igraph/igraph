
context("multilevel.community")

test_that("multilevel.community works", {

  library(igraph)

  g <- graph.famous("Zachary")
  mc <- multilevel.community(g)
  
  expect_that(membership(mc),
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
