
context("transitivity")

test_that("transitivity works", {
  library(igraph)
  set.seed(42)
  g <- erdos.renyi.game(100, p=10/100)

  t1 <- transitivity(g, type="global")
  expect_that(t1, equals(0.10483870967741935887))

  t2 <- transitivity(g, type="average")
  expect_that(t2, equals(0.10159943848720931481))

  t3 <- transitivity(g, type="local", vids=V(g))
  t33 <- transitivity(g, type="local")
  est3 <- structure(c(0, 0.06667, 0.1028, 0.1016, 0.1333, 0.2222),
                    .Names = c("Min.", "1st Qu.", "Median", "Mean",
                      "3rd Qu.", "Max."),
                    class = c("summaryDefault", "table"))
  expect_that(summary(t3), equals(est3))
  expect_that(summary(t33), equals(est3))
})
