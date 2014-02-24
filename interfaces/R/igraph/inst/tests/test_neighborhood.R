
context("neighborhood")

test_that("neighborhood works", {

  library(igraph)

  neig <- function(graph, order, vertices) {
    sp <- shortest.paths(graph)
    v <- unique(unlist(lapply(vertices, function(x) {
      w <- which(sp[x,] <= order)
    })))
    induced.subgraph(graph, c(v,vertices))
  }

  g <- erdos.renyi.game(50, 5/50)

  v <- sample(vcount(g), 1)
  g1 <- graph.neighborhood(g, 2, v)[[1]]
  g2 <- neig(g, 2, v)
  expect_that(graph.isomorphic(g1, g2), is_true())

#########

  nei <- function(graph, order, vertices) {
    sp <- shortest.paths(graph)
    v <- unique(unlist(lapply(vertices, function(x) {
      w <- which(sp[x,] <= order)
    })))
    v
  }

  v1 <- neighborhood(g, 2, v)[[1]]
  v2 <- nei(g, 2, v)
  expect_that(sort(v1), equals(sort(v2))) 

#########

  s <- neighborhood.size(g, 2, v)[[1]]
  expect_that(s, equals(length(v1)))

})

test_that("mindist works", {

  library(igraph)
  g <- graph.ring(10)
  expect_that(neighborhood.size(g, order=2, mindist=0), equals(rep(5, 10)))
  expect_that(neighborhood.size(g, order=2, mindist=1), equals(rep(4, 10)))
  expect_that(neighborhood.size(g, order=2, mindist=2), equals(rep(2, 10)))

  n0 <- neighborhood(g, order=2, 5:6, mindist=0)
  n1 <- neighborhood(g, order=2, 5:6, mindist=1)
  n2 <- neighborhood(g, order=2, 5:6, mindist=2)

  expect_that(lapply(n0, sort), equals(list(3:7, 4:8)))
  expect_that(lapply(n1, sort), equals(list(c(3,4,6,7), c(4,5,7,8))))
  expect_that(lapply(n2, sort), equals(list(c(3,7), c(4,8))))

  ng0 <- graph.neighborhood(g, order=2, 5:6, mindist=0)
  ng1 <- graph.neighborhood(g, order=2, 5:6, mindist=1)
  ng2 <- graph.neighborhood(g, order=2, 5:6, mindist=2)

  expect_that(sapply(ng0, vcount), equals(c(5,5)))
  expect_that(sapply(ng1, vcount), equals(c(4,4)))
  expect_that(sapply(ng2, vcount), equals(c(2,2)))

  expect_that(sapply(ng0, ecount), equals(c(4,4)))
  expect_that(sapply(ng1, ecount), equals(c(2,2)))
  expect_that(sapply(ng2, ecount), equals(c(0,0)))

})
