
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

test_that("weighted optimal.community works", {

  library(igraph)
  set.seed(42)
  g <- graph.full(5) + graph.ring(5)
  E(g)$weight <- sample(1:2, ecount(g), replace=TRUE)

  oc <- optimal.community(g)
  expect_that(modularity(oc), equals(0.4032))
})
