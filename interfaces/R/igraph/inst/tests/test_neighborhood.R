
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
