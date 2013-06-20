
context("infix operators")

test_that("infix operators work", {

  library(igraph)

  g <- graph.ring(10)
  V(g)$name <- letters[1:10]
  E(g)$name <- LETTERS[1:10]

  g <- g - c("a", "b")
  expect_that(vcount(g), equals(8))
  expect_that(ecount(g), equals(7))
  expect_that(graph.isomorphic(g, graph.lattice(8)), is_true())

  g <- g - edge("e|f")
  expect_that(graph.isomorphic(g, graph.lattice(5) + graph.lattice(3)),
              is_true())

  g <- g - edge("H")
  expect_that(graph.isomorphic(g, graph.formula(a-b-c, d-e-f, g-h)),
              is_true())

  g <- graph.ring(10)
  V(g)$name <- letters[1:10]

  g <- g - path("a", "b", "c", "d")
  expect_that(graph.isomorphic(g, graph.lattice(8) + 2), is_true())

  expect_that(graph.isomorphic(g - V(g)[c('d', 'g')],
                               graph.lattice(4) + graph.lattice(2) + 2),
              is_true())

  expect_that(graph.isomorphic(g - E(g)['f' %--% 'g'],
                               graph.lattice(5) + graph.lattice(3) + 2),
              is_true())

})
