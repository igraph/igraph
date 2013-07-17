
context("farthest.nodes")

test_that("farthest.nodes works", {

  library(igraph)

  kite <- graph.formula(Andre    - Beverly:Carol:Diane:Fernando,
                        Beverly  - Andre:Diane:Ed:Garth,
                        Carol    - Andre:Diane:Fernando,
                        Diane    - Andre:Beverly:Carol:Ed:Fernando:Garth,
                        Ed       - Beverly:Diane:Garth,
                        Fernando - Andre:Carol:Diane:Garth:Heather,
                        Garth    - Beverly:Diane:Ed:Fernando:Heather,
                        Heather  - Fernando:Garth:Ike,
                        Ike      - Heather:Jane,
                        Jane     - Ike)

  fn <- farthest.nodes(kite)
  expect_that(fn, equals(c(1,10,4)))

  expect_that(shortest.paths(kite, v=fn[1], to=fn[2])[1], equals(fn[3]))
  expect_that(diameter(kite), equals(fn[3]))

})
