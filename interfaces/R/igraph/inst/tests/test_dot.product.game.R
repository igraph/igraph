
context("Dot-product random graphs")

test_that("Dot product rng works", {

  library(igraph)
  set.seed(42)
  vecs <- cbind(c(0,1,1,1,0)/3, c(0,1,1,0,1)/3, c(1,1,1,1,0)/4,
                c(0,1,1,1,0))

  g <- dot.product.game(vecs)
  g0 <- graph.formula(1:2:3-4)
  expect_that(g[], is_equivalent_to(g0[]))

  g2 <- dot.product.game(vecs, directed=TRUE)
  g20 <- graph.formula(1:2:3:4, 1-+3, 1-+4, 3-+4, 4+-1, 4+-3)
  expect_that(g[], is_equivalent_to(g20[]))

  vecs <- replicate(5, rep(1/2, 4))
  g <- dot.product.game(vecs)
  expect_that(g[], is_equivalent_to(graph.full(5)[]))

  g2 <- dot.product.game(vecs, directed=TRUE)
  expect_that(g2[], is_equivalent_to(graph.full(5, directed=TRUE)[]))

  vecs <- replicate(100, rep(sqrt(1/8), 4))
  g <- dot.product.game(vecs)
  expect_that(ecount(g), equals(2454))

  g2 <- dot.product.game(vecs, directed=TRUE)
  expect_that(ecount(g2), equals(4938))
  
})
