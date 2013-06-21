
context("graph.knn")

test_that("graph.knn works", {
  library(igraph)
  set.seed(42)

  ## Some trivial ones
  g <- graph.ring(10)
  expect_that(graph.knn(g), equals(list(knn=rep(2,10), knnk=c(NaN, 2))))

  g2 <- graph.star(10)
  expect_that(graph.knn(g2), equals(list(knn=c(1, rep(9,9)),
                                         knnk=c(9, rep(NaN, 7), 1))))

  ## A scale-free one, try to plot 'knnk'
  g3 <- simplify(ba.game(1000, m=5))
  r3 <- graph.knn(g3)
  expect_that(r3$knn[43], equals(46))
  expect_that(r3$knn[1000], equals(192.4))
  expect_that(r3$knnk[100], equals(18.78))
  expect_that(length(r3$knnk), equals(359))
  
  ## A random graph
  g4 <- random.graph.game(1000, p=5/1000)
  r4 <- graph.knn(g4)
  expect_that(r4$knn[1000], equals(20/3))
  expect_that(length(r4$knnk), equals(15))
  expect_that(r4$knnk[12], equals(19/3))

  ## A weighted graph
  g5 <- graph.star(10)
  E(g5)$weight <- seq(ecount(g5))
  r5 <- graph.knn(g5)
  expect_that(r5, equals(structure(list(knn = c(1, 45, 22.5, 15,
                                          11.25, 9, 7.5, 6.42857142857143,
                                          5.625, 5), knnk =
                                        c(14.1448412698413, NaN, NaN,
                                          NaN, NaN, NaN, NaN, NaN,
                                          1)), .Names = c("knn",
                                                 "knnk")) ))
})
