
context("add_vertices")

test_that("add_vertices works", {
  library(igraph)
  g <- graph_from_literal(A-B-C-D-E)
  g2 <- add_vertices(g, (nv <- 4))
  expect_that(vcount(g2), equals(vcount(g) + nv))
  expect_that(ecount(g2), equals(ecount(g)))
  expect_that(as_edgelist(g2), equals(as_edgelist(g)))
})

test_that("add_vertices handles attributes properly", {
  library(igraph)
  g <- graph_from_literal(A-B-C-D-E)
  g3 <- add_vertices(g, (nv <- 3),
                     attr=list(name=(names <- c("F","G","H")),
                       weight=weights <- 1:3))
  expect_that(V(g3)$name, equals(c(V(g)$name, names)))
  expect_that(V(g3)$weight, equals(c(rep(NA, vcount(g)), weights)))
})
