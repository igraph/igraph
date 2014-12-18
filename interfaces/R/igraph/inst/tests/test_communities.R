
context("communities")

test_that("community detection functions work", {
  library(igraph)
  set.seed(42)

  F <- list("cluster_edge_betweenness", "cluster_fast_greedy",
            "cluster_label_prop", "cluster_leading_eigen",
            "cluster_louvain", "cluster_optimal",
            "cluster_spinglass", "cluster_walktrap")

  karate <- make_graph("Zachary")

  for (f in F) {
    f <- get(f)
    comm <- f(karate)
    
    expect_that(modularity(comm),
                equals(modularity(karate, membership(comm))))
  
    cc <- communities(comm)
    expect_that(all(!duplicated(unlist(cc))), is_true())
    expect_that(all(unlist(cc) <= vcount(karate) & unlist(cc) >= 1),
                is_true())
    expect_that(length(comm), equals(max(membership(comm))))
  }

  fc <- cluster_fast_greedy(karate)
  m1 <- modularity(karate, cut_at(fc, no=1))
  m2 <- modularity(karate, cut_at(fc, no=2))
  m3 <- modularity(karate, cut_at(fc, no=3))
  m4 <- modularity(karate, cut_at(fc, no=4))
  expect_that(m1, equals(0))
  expect_that(m2, equals(0.3717948718))
  expect_that(m3, equals(0.3806706114))
  expect_that(m4, equals(0.3759861933))

  cr <- crossing(fc, karate)
  expect_that(cr, equals(c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE,
                           TRUE, FALSE, FALSE, TRUE, TRUE, TRUE,
                           FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
                           FALSE, TRUE, FALSE, TRUE, FALSE, FALSE,
                           TRUE, TRUE, TRUE, FALSE, TRUE, FALSE,
                           FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                           FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                           TRUE, TRUE, FALSE, FALSE, FALSE, FALSE,
                           FALSE, FALSE, TRUE, FALSE, FALSE, FALSE,
                           FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                           FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                           FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                           FALSE, FALSE, FALSE, FALSE) ))
})

test_that("creating communities objects works", {
  library(igraph)
  set.seed(42)

  karate <- make_graph("Zachary")

  membership <- sample(1:2, vcount(karate), replace=TRUE)
  mod <- modularity(karate, membership)
  comm <- make_clusters(algorithm="random", membership=membership,
                             modularity = mod)

  expect_that(as.vector(membership(comm)), equals(membership))
  expect_that(modularity(comm), equals(mod))
  expect_that(algorithm(comm), equals("random"))

})

test_that("communities function works", {
  library(igraph)
  g <- make_graph("Zachary")
  oc <- cluster_optimal(g)
  gr <- communities(oc)
  expect_that(gr, equals
    (structure(list(`1` = c(1L, 2L, 3L, 4L, 8L, 12L, 13L, 14L, 18L,
     20L, 22L), `2` = c(5L, 6L, 7L, 11L, 17L), `3` = c(9L, 10L, 15L,
     16L, 19L, 21L, 23L, 27L, 30L, 31L, 33L, 34L), `4` = c(24L, 25L,
     26L, 28L, 29L, 32L)), .Dim = 4L, .Dimnames = list(c("1", "2",
     "3", "4")))))

  g <- make_ring(5) + make_ring(5)
  V(g)$name <- letters[1:10]
  oc <- cluster_optimal(g)
  gr <- communities(oc)
  expect_that(gr, equals(structure(list(`1` = letters[1:5],
                                        `2` = letters[6:10]),
                                        .Dim = 2L,
                                        .Dimnames = list(c("1", "2")))))
})
