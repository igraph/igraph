
context("Local scan statistics")

library(igraph)
require(digest)

set.seed(12345)
n <- 10^3
p <- 0.1
g <- erdos.renyi.game(n,p)
E(g)$weight = sample(ecount(g))
gp <- erdos.renyi.game(n,p)
E(gp)$weight = sample(ecount(gp))

test_that("General scan-stat works, US, scan-0, unweighted", {
  s1 <- local.scan(g, k=0)
  expect_that(digest(s1), equals("659ffaaf303742f0806a79b8ff3d88b3"))
})

test_that("General scan-stat works, US, scan-0, weighted", {
  s1 <- local.scan(g, k=0, weighted=TRUE)
  expect_that(digest(s1), equals("0f8d7ac831389cea04e0bfc5e2510c73"))
})


test_that("General scan-stat works, US, scan-1, unweighted", {
  s1 <- local.scan(g)
  expect_that(digest(s1), equals("df0fd77489f70cc47f682dc31d9f52f5"))
})

test_that("General scan-stat works, US, scan-1, weighted", {
  s1 <- local.scan(g, k=1, weighted=TRUE)
  expect_that(digest(s1), equals("af720916ae4b49881745d2dcdd614401"))
})

test_that("General scan-stat works, US, scan-2, unweighted", {
  s1 <- local.scan(g, k=2)
  expect_that(digest(s1), equals("6f47f47abde25d00d615dd56826cca5a"))
})

test_that("General scan-stat works, US, scan-2, weighted", {
  s1 <- local.scan(g, k=2, weighted=TRUE)
  expect_that(digest(s1), equals("e02e9d58168ee5d53850497f6d4c76b0"))
})

test_that("General scan-stat works, THEM, scan-0, unweighted", {
  s1 <- local.scan(g, gp, k=0)
  expect_that(digest(s1), equals("f584f7d287f8f89f5f7882165ca41b8c"))
})

test_that("General scan-stat works, THEM, scan-0, weighted", {
  s1 <- local.scan(g, gp, k=0, weighted=TRUE)
  expect_that(digest(s1), equals("213db8e7517d1e6406da3dbd55281ed1"))
})

test_that("General scan-stat works, THEM, scan-1, unweighted", {
  s1 <- local.scan(g, gp, k=1)
  expect_that(digest(s1), equals("e9ca740ebba2fd1db4abe939954b2638"))
})

test_that("General scan-stat works, THEM, scan-1, weighted", {
  s1 <- local.scan(g, gp, k=1, weighted=TRUE)
  expect_that(digest(s1), equals("a98e9a03eda7feaae8524dc9348ad74b"))
})

test_that("General scan-stat works, THEM, scan-2, unweighted", {
  s1 <- local.scan(g, gp, k=2)
  expect_that(digest(s1), equals("a3237a9a55e9d86ab471c81a291eb03b"))
})

test_that("General scan-stat works, THEM, scan-2, weighted", {
  s1 <- local.scan(g, gp, k=2, weighted=TRUE)
  expect_that(digest(s1), equals("995d0b6a952834ff6e534efc2cfb917b"))
})


set.seed(42)
n <- 10^3
p <- 0.1
g <- erdos.renyi.game(n, p, directed=TRUE)
E(g)$weight = sample(ecount(g))
gp <- erdos.renyi.game(n, p)
E(gp)$weight = sample(ecount(gp))

## US, scan-0, unweighted, directed
## TODO

test_that("General scan-stat works, US, scan-1, unweighted, directed", {

  s1o <- local.scan(g, k=1, weighted=FALSE, mode="out")
  expect_that(digest(s1o), equals("ac463c21b2b6bc91abf82f0141a4a7d4"))

  s1i <- local.scan(g, k=1, weighted=FALSE, mode="in")
  expect_that(digest(s1i), equals("13fdaaeec54118e217821b56d8c3ff03"))

})

test_that("General scan-stat works, US, scan-1, weighted, directed", {

  s1o <- local.scan(g, k=1, weighted=TRUE, mode="out")
  expect_that(digest(s1o), equals("da8e14f2ba63efc74b5fd7b9d8f79bbc"))

  s1i <- local.scan(g, k=1, weighted=TRUE, mode="in")
  expect_that(digest(s1i), equals("f5f07eebb907ae0a244195a20971be11"))

})

## US, scan-2, unweighted, directed
## TODO

test_that("Local scan-1 approximations work", {
  
  set.seed(42)

  for (i in 1:100) {

    g <- erdos.renyi.game(100, 10/100)
    
    s1 <- local.scan(g)
    s1a <- local.scan1.ecount.approx(g, 20)$res
    expect_that(cor(s1, s1a) > 0.95, is_true())

    E <- graph.eigen(g, which=list(howmany=20, pos="LM"))
    s1aa <- local.scan1.ecount.approx.eigen(g, E$values, E$vectors)
    expect_that(cor(s1, s1aa) > 0.95, is_true())

    E2 <- eigen(get.adjacency(g, sparse=FALSE))
    s1aaa <- colSums(E2$values ^3 * t(E2$vectors)^2 / 2) + degree(g)
    expect_that(s1aaa, equals(s1))

  }

})

test_that("Issue 18 is resolved", {

  library(igraph)
  g <- graph(c(1,2,2,1, 1,3,3,1, 2,4, 3,4, 3,5,5,3, 4,5,5,4))
  expect_that(local.scan(g, mode="all"), equals(c(4, 3, 7, 6, 5)))
  expect_that(local.scan(g, mode="out"), equals(c(4, 3, 7, 2, 5)))
  expect_that(local.scan(g, mode="in"), equals(c(4, 2, 4, 6, 5)))
})
