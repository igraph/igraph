
context("Scan-1")

test_that("scan-1 statistics works", {
  
  library(igraph)
  set.seed(42)

  for (i in 1:100) {

    g <- erdos.renyi.game(100, 10/100)
    
    s1 <- scan1(g)
    s1a <- scan1.approx(g, 20)$res
    expect_that(cor(s1, s1a) > 0.95, is_true())

    E <- graph.eigen(g, which=list(howmany=20, pos="LM"))
    s1aa <- scan1.approx.eigen(g, E$values, E$vectors)
    expect_that(cor(s1, s1aa) > 0.95, is_true())

    E2 <- eigen(get.adjacency(g, sparse=FALSE))
    s1aaa <- colSums(E2$values ^3 * t(E2$vectors)^2 / 2) + degree(g)
    expect_that(s1aaa, equals(s1))

  }

})
