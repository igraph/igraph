
context("Hierarchical random graphs")

test_that("Starting from state works (#225)", {
  library(igraph)
  set.seed(42)
  
  g <- erdos.renyi.game(10, p=1/2) + erdos.renyi.game(10, p=1/2)
  hrg <- hrg.fit(g)
  hrg2 <- hrg.fit(g, hrg=hrg, start=TRUE, steps=1)
  expect_that(hrg2, is_equivalent_to(hrg))
  
})
