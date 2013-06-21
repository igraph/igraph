
context("optimal.community")

test_that("optimal.community works", {

  library(igraph)
  g <- graph.famous("Zachary")
  oc <- optimal.community(g)

  expect_that(membership(oc),
              equals(c(1, 1, 1, 1, 2, 2, 2, 1, 3, 3, 2, 1, 1, 1, 3, 3,
                       2, 1, 3, 1, 3, 1, 3, 4, 4, 4, 3, 4, 4, 3, 3, 4,
                       3, 3) ))
  expect_that(modularity(g, oc$membership), equals(oc$modularity))
  expect_that(length(oc), equals(4))
  expect_that(sizes(oc),
              equals(structure(c(11L, 5L, 12L, 6L), .Dim=4L,
              .Dimnames=structure(list(`Community sizes`=c("1", "2",
              "3", "4")), .Names="Community sizes"), class="table") ))

})
