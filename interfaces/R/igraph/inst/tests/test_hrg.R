
context("Hierarchical random graphs")

test_that("Starting from state works (#225)", {
  library(igraph)
  set.seed(42)

  res <- structure(list(left = c(-12, 13, -19, 4, 15, -3, 0, -13, -15, 
                          -18, -6, -14, -7, -9, -8, 2, -10, 10, 11),
                        right = c(-17, -5, 17, 7, 16, -2, 1, -4, 5,
                          14, 12, 8, 3, 6, -16, 9, 18, -11, 19),
                        prob = c(0, 0.5, 0, 0, 0, 0.888888888888889,
                          0, 0.5, 0.714285714285714, 0.375,
                          0.333333333333333, 0.444444444444444, 1,
                          0.625, 0.8, 0,  0.222222222222222,
                          0.714285714285714, 1),
                        edges = c(0, 1, 0, 0, 0, 8, 0, 3, 5, 3, 2, 4,
                          2, 5, 8, 0, 2, 5, 1),
                        vertices = c(20, 3, 3, 2, 2, 6, 2, 5, 8, 9, 7,
                          10, 3, 9, 7, 2, 10, 8, 2)),
                   .Names = c("left", "right", "prob", "edges",
                     "vertices"),
                   class = "igraphHRG")
  
  g <- erdos.renyi.game(10, p=1/2) + erdos.renyi.game(10, p=1/2)
  hrg <- hrg.fit(g)
  hrg

  expect_that(hrg, is_equivalent_to(res))
  
  hrg2 <- hrg.fit(g, hrg=hrg, start=TRUE, steps=1)
  expect_that(hrg2, is_equivalent_to(res))
  
})
