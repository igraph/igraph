
context("degree")

test_that("degree works", {
  library(igraph)
  
  g <- erdos.renyi.game(100, 1/100)
  d <- degree(g)
  el <- get.edgelist(g)
  expect_that(as.numeric(table(el)), equals(d[d!=0]))

  expect_that(degree(g) / (vcount(g)-1), equals(degree(g, normalized=TRUE)))

  g2 <- erdos.renyi.game(100, 2/100, dir=TRUE)
  din <- degree(g2, mode="in")
  dout <- degree(g2, mode="out")
  el2 <- get.edgelist(g2)
  expect_that(as.numeric(table(el2[,1])), equals(dout[dout!=0]))
  expect_that(as.numeric(table(el2[,2])), equals(din[din!=0]))

  expect_that(degree(g2, mode="in") / (vcount(g2)-1),
              equals(degree(g2, mode="in", normalized=TRUE)))
  expect_that(degree(g2, mode="out") / (vcount(g2)-1),
              equals(degree(g2, mode="out", normalized=TRUE)))
  expect_that(degree(g2, mode="all") / (vcount(g2)-1),
              equals(degree(g2, mode="all", normalized=TRUE)))
})
