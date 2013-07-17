
context("Bug 1073705")

test_that("Weighted indexing does not remove edges", {
  library(igraph)

  g <- graph.ring(10)
  g[1, 2, attr="weight"] <- 0
  expect_that("weight" %in% list.edge.attributes(g), is_true())
  expect_that(E(g)$weight, equals(c(0, rep(NA, 9))))

  el <- get.edgelist(g)
  g[from=el[,1], to=el[,2], attr="sim"] <- rep(0:1, length=ecount(g))
  expect_that("sim" %in% list.edge.attributes(g), is_true())
  expect_that(E(g)$sim, equals(rep(0:1, 5)))

  V(g)$name <- letters[seq_len(vcount(g))]
  el <- get.edgelist(g)
  g[from=el[,1], to=el[,2], attr="sim"] <- rep(1:0, length=ecount(g))
  expect_that(E(g)$sim, equals(rep(1:0, 5)))
})


