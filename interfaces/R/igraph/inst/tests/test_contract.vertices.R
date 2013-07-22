
context("contract.vertices")

test_that("contract.vertices works", {
  library(igraph)
  set.seed(42)

  g <- graph.ring(10)
  g$name <- "Ring"
  V(g)$name <- letters[1:vcount(g)]
  E(g)$weight <- sample(ecount(g))

  g2 <- contract.vertices(g, rep(1:5, each=2),
                          vertex.attr.comb=toString)

  ## graph and edge attributes are kept, vertex attributes are
  ## combined using the 'toString' function.
  expect_that(g2$name, equals(g$name))
  expect_that(V(g2)$name, equals(c("a, b", "c, d", "e, f", "g, h", "i, j")))
  expect_that(as.matrix(g2[]),
              is_equivalent_to(cbind(c(10,9,0,0,7), c(9,3,6,0,0),
                                     c(0,6,4,8,0), c(0,0,8,5,1),
                                     c(7,0,0,1,2))))
})
