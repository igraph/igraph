
context("Graphlets")

test_that("Getting subcliques works", {
  library(igraph)
  set.seed(42*42)
  g <- erdos.renyi.game(10, 4/10)
  E(g)$weight <- as.double(sample(1:10, ecount(g), replace=TRUE))
  ids <- 1:vcount(g)

  cl <- maximal.cliques(g)
  cl <- lapply(cl, "-", 1)[c(9, 2, 3, 10, 5, 7, 6, 1, 4, 8)]

  res <- .Call("R_igraph_subclique_next", g, E(g)$weight, ids, cl,
               PACKAGE="igraph")

  for (i in seq_along(res$graphs)) {
    V(res$graphs[[i]])$name <- res$ids[[i]]
    E(res$graphs[[i]])$weight <- res$weights[[i]]    
  }
  
  expect_that(res$thr, equals(c(7,2,5,3,1,7,3,2,4,7)))
  expect_that(res$next_thr, equals(c(Inf, 4, 8, 8, Inf, 9, 5, 3, 5, Inf)))
  expect_that(res$weights, equals(list(numeric(), c(4,8), c(8,8), c(9,8),
                                       numeric(), c(9,9), c(7,7,9,10,5),
                                       c(7,3,4,5), c(5,7,5,9,10),
                                       numeric())))
  expect_that(res$ids, equals(list(integer(), c(5,9,10), c(1,9,10),
                                   c(1,8,9), integer(), c(3,7,6),
                                   c(3,6,5,4), c(3,5,10,2), c(2,5,3,4),
                                   integer())))
  expect_that(sapply(res$graphs, vcount), equals(sapply(res$ids, length)))
  expect_that(sapply(res$graphs, ecount),
              equals(sapply(res$weights, length)))
})
