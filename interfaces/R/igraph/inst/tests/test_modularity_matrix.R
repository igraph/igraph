
context("mod.matrix")

test_that("mod.matrix works", {

  library(igraph)

  kar <- graph.famous("zachary")

  fc <- fastgreedy.community(kar)

  m1 <- modularity(kar, membership(fc))
  m2 <- modularity(kar, membership(fc), weights=rep(1, ecount(kar)))
  expect_that(m1, equals(m2))
  
  B1 <- mod.matrix(kar, membership(fc))
  B2 <- mod.matrix(kar, membership(fc), weights=rep(1, ecount(kar)))

  expect_that(B1, equals(B2))

})
