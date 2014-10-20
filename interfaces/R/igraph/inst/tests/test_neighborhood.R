
context("ego")

test_that("ego works", {

  library(igraph)

  neig <- function(graph, order, vertices) {
    sp <- distances(graph)
    v <- unique(unlist(lapply(vertices, function(x) {
      w <- which(sp[x,] <= order)
    })))
    induced_subgraph(graph, c(v,vertices))
  }

  g <- sample_gnp(50, 5/50)

  v <- sample(vcount(g), 1)
  g1 <- make_ego_graph(g, 2, v)[[1]]
  g2 <- neig(g, 2, v)
  expect_that(graph.isomorphic(g1, g2), is_true())

#########

  nei <- function(graph, order, vertices) {
    sp <- distances(graph)
    v <- unique(unlist(lapply(vertices, function(x) {
      w <- which(sp[x,] <= order)
    })))
    v
  }

  v1 <- ego(g, 2, v)[[1]]
  v2 <- nei(g, 2, v)
  expect_that(as.vector(sort(v1)), equals(sort(v2)))

#########

  s <- ego_size(g, 2, v)[[1]]
  expect_that(s, equals(length(v1)))

})

test_that("mindist works", {

  library(igraph)
  g <- make_ring(10)
  expect_that(ego_size(g, order=2, mindist=0), equals(rep(5, 10)))
  expect_that(ego_size(g, order=2, mindist=1), equals(rep(4, 10)))
  expect_that(ego_size(g, order=2, mindist=2), equals(rep(2, 10)))

  unvs <- function(x) lapply(x, as.vector)

  n0 <- unvs(ego(g, order=2, 5:6, mindist=0))
  n1 <- unvs(ego(g, order=2, 5:6, mindist=1))
  n2 <- unvs(ego(g, order=2, 5:6, mindist=2))

  expect_that(lapply(n0, sort), equals(list(3:7, 4:8)))
  expect_that(lapply(n1, sort), equals(list(c(3,4,6,7), c(4,5,7,8))))
  expect_that(lapply(n2, sort), equals(list(c(3,7), c(4,8))))

  ng0 <- make_ego_graph(g, order=2, 5:6, mindist=0)
  ng1 <- make_ego_graph(g, order=2, 5:6, mindist=1)
  ng2 <- make_ego_graph(g, order=2, 5:6, mindist=2)

  expect_that(sapply(ng0, vcount), equals(c(5,5)))
  expect_that(sapply(ng1, vcount), equals(c(4,4)))
  expect_that(sapply(ng2, vcount), equals(c(2,2)))

  expect_that(sapply(ng0, ecount), equals(c(4,4)))
  expect_that(sapply(ng1, ecount), equals(c(2,2)))
  expect_that(sapply(ng2, ecount), equals(c(0,0)))

})
