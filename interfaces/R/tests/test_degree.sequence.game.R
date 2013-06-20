
context("degree.sequence.game")

test_that("degree.sequence.game works", {
  library(igraph)

  gc <- function(graph) {
    clu <- clusters(graph)
    induced.subgraph(graph, which(clu$membership==which.max(clu$csize)))
  }

  g <- gc(erdos.renyi.game(1000, 2/1000))

  nG <- degree.sequence.game(degree(g), method="simple")
  expect_that(degree(nG), equals(degree(g)))

  nG <- degree.sequence.game(degree(g), method="vl")
  expect_that(degree(nG), equals(degree(g)))
  expect_that(is.connected(nG), is_true())
  expect_that(is.simple(nG), is_true())

  #####

  g <- erdos.renyi.game(1000, 1/1000)

  nG <- degree.sequence.game(degree(g), method="simple")
  expect_that(degree(nG), equals(degree(g)))

  g2 <- erdos.renyi.game(1000, 2/1000, dir=TRUE)

  nG2 <- degree.sequence.game(degree(g, mode="out"),
                              degree(g, mode="in"),
                              method="simple")
  expect_that(degree(nG, mode="out"), equals(degree(g, mode="out")))
  expect_that(degree(nG, mode="in"), equals(degree(g, mode="in")))
})
