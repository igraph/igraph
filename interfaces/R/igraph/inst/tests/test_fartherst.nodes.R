
context("farthest_vertices")

test_that("farthest_vertices works", {

  library(igraph)

  kite <- graph_from_literal(Andre    - Beverly:Carol:Diane:Fernando,
                    Beverly  - Andre:Diane:Ed:Garth,
                    Carol    - Andre:Diane:Fernando,
                    Diane    - Andre:Beverly:Carol:Ed:Fernando:Garth,
                    Ed       - Beverly:Diane:Garth,
                    Fernando - Andre:Carol:Diane:Garth:Heather,
                    Garth    - Beverly:Diane:Ed:Fernando:Heather,
                    Heather  - Fernando:Garth:Ike,
                    Ike      - Heather:Jane,
                    Jane     - Ike)

  fn <- farthest_vertices(kite)
  fn$vertices <- as.vector(fn$vertices)
  expect_that(fn, equals(list(vertices = c(1, 10), distance = 4)))

  expect_that(distances(kite, v=fn$vertices[1], to=fn$vertices[2])[1],
              equals(fn$distance))
  expect_that(diameter(kite), equals(fn$distance))

})
