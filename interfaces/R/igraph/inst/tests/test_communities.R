
context("communities")

test_that("community detection functions work", {
  library(igraph)
  set.seed(42)

  F <- list("edge.betweenness.community", "fastgreedy.community",
            "label.propagation.community", "leading.eigenvector.community",
            "multilevel.community", "optimal.community",
            "spinglass.community", "walktrap.community")

  karate <- graph.famous("Zachary")

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

  fc <- fastgreedy.community(karate)
  m1 <- modularity(karate, cutat(fc, no=1))
  m2 <- modularity(karate, cutat(fc, no=2))
  m3 <- modularity(karate, cutat(fc, no=3))
  m4 <- modularity(karate, cutat(fc, no=4))
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

  karate <- graph.famous("Zachary")

  membership <- sample(1:2, vcount(karate), replace=TRUE)
  mod <- modularity(karate, membership)
  comm <- create.communities(algorithm="random", membership=membership,
                             mod=mod, foo="bar")
  print(comm)

  expect_that(membership(comm), equals(membership))
  expect_that(modularity(comm), equals(mod))
  expect_that(algorithm(comm), equals("random"))
  expect_that(comm$foo, equals("bar"))

})
