
context("cliques")

test_that("cliques works", {
  library(igraph)
  set.seed(42)

  check.clique <- function(graph, vids) {
    s <- induced.subgraph(graph, vids)
    ecount(s) == vcount(s) * (vcount(s)-1) / 2
  }

  g <- erdos.renyi.game(100, 0.3)
  expect_that(clique.number(g), equals(6))
  
  cl <- sapply(cliques(g, min=6), check.clique, graph=g)
  lcl <- sapply(largest.cliques(g), check.clique, graph=g)
  expect_that(cl, equals(lcl))
  expect_that(cl, equals(rep(TRUE, 17)))
  expect_that(lcl, equals(rep(TRUE, 17)))

  ## To have a bit less maximal cliques, about 100-200 usually
  g <- erdos.renyi.game(100, 0.03)
  expect_that(all(sapply(maximal.cliques(g), check.clique, graph=g)),
              is_true())
})
