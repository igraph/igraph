
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
  s1 <- local_scan(g, k=0)
  expect_that(digest(s1), equals("659ffaaf303742f0806a79b8ff3d88b3"))
})

test_that("General scan-stat works, US, scan-0, weighted", {
  s1 <- local_scan(g, k=0, weighted=TRUE)
  expect_that(digest(s1), equals("0f8d7ac831389cea04e0bfc5e2510c73"))
})


test_that("General scan-stat works, US, scan-1, unweighted", {
  s1 <- local_scan(g)
  expect_that(digest(s1), equals("df0fd77489f70cc47f682dc31d9f52f5"))
})

test_that("General scan-stat works, US, scan-1, weighted", {
  s1 <- local_scan(g, k=1, weighted=TRUE)
  expect_that(digest(s1), equals("af720916ae4b49881745d2dcdd614401"))
})

test_that("General scan-stat works, US, scan-2, unweighted", {
  s1 <- local_scan(g, k=2)
  expect_that(digest(s1), equals("6f47f47abde25d00d615dd56826cca5a"))
})

test_that("General scan-stat works, US, scan-2, weighted", {
  s1 <- local_scan(g, k=2, weighted=TRUE)
  expect_that(digest(s1), equals("e02e9d58168ee5d53850497f6d4c76b0"))
})

test_that("General scan-stat works, THEM, scan-0, unweighted", {
  s1 <- local_scan(g, gp, k=0)
  expect_that(digest(s1), equals("f584f7d287f8f89f5f7882165ca41b8c"))
})

test_that("General scan-stat works, THEM, scan-0, weighted", {
  s1 <- local_scan(g, gp, k=0, weighted=TRUE)
  expect_that(digest(s1), equals("213db8e7517d1e6406da3dbd55281ed1"))
})

test_that("General scan-stat works, THEM, scan-1, unweighted", {
  s1 <- local_scan(g, gp, k=1)
  expect_that(digest(s1), equals("e9ca740ebba2fd1db4abe939954b2638"))
})

test_that("General scan-stat works, THEM, scan-1, weighted", {
  s1 <- local_scan(g, gp, k=1, weighted=TRUE)
  expect_that(digest(s1), equals("a98e9a03eda7feaae8524dc9348ad74b"))
})

test_that("General scan-stat works, THEM, scan-2, unweighted", {
  s1 <- local_scan(g, gp, k=2)
  expect_that(digest(s1), equals("a3237a9a55e9d86ab471c81a291eb03b"))
})

test_that("General scan-stat works, THEM, scan-2, weighted", {
  s1 <- local_scan(g, gp, k=2, weighted=TRUE)
  expect_that(digest(s1), equals("995d0b6a952834ff6e534efc2cfb917b"))
})

test_that("Neighborhoods work for us", {
  nei <- neighborhood(g, order=1)
  s1 <- local_scan(g, neighborhoods=nei)
  expect_that(digest(s1), equals("df0fd77489f70cc47f682dc31d9f52f5"))
  s1 <- local_scan(g, k=1, weighted=TRUE, neighborhoods=nei)
  expect_that(digest(s1), equals("af720916ae4b49881745d2dcdd614401"))

  nei <- neighborhood(g, order=2)
  s1 <- local_scan(g, k=2, neighborhoods=nei)
  expect_that(digest(s1), equals("6f47f47abde25d00d615dd56826cca5a"))
  s1 <- local_scan(g, k=2, weighted=TRUE, neighborhoods=nei)
  expect_that(digest(s1), equals("e02e9d58168ee5d53850497f6d4c76b0"))

})

test_that("Neighborhoods work for them", {

  nei <- neighborhood(g, order=1)
  s1 <- local_scan(g, gp, k=1, neighborhoods=nei)
  expect_that(digest(s1), equals("e9ca740ebba2fd1db4abe939954b2638"))
  s1 <- local_scan(g, gp, k=1, weighted=TRUE, neighborhoods=nei)
  expect_that(digest(s1), equals("a98e9a03eda7feaae8524dc9348ad74b"))

  nei <- neighborhood(g, order=2)
  s1 <- local_scan(g, gp, k=2, neighborhoods=nei)
  expect_that(digest(s1), equals("a3237a9a55e9d86ab471c81a291eb03b"))
  s1 <- local_scan(g, gp, k=2, weighted=TRUE, neighborhoods=nei)
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

  s1o <- local_scan(g, k=1, weighted=FALSE, mode="out")
  expect_that(digest(s1o), equals("ac463c21b2b6bc91abf82f0141a4a7d4"))

  s1i <- local_scan(g, k=1, weighted=FALSE, mode="in")
  expect_that(digest(s1i), equals("13fdaaeec54118e217821b56d8c3ff03"))

})

test_that("General scan-stat works, US, scan-1, weighted, directed", {

  s1o <- local_scan(g, k=1, weighted=TRUE, mode="out")
  expect_that(digest(s1o), equals("da8e14f2ba63efc74b5fd7b9d8f79bbc"))

  s1i <- local_scan(g, k=1, weighted=TRUE, mode="in")
  expect_that(digest(s1i), equals("f5f07eebb907ae0a244195a20971be11"))

})

## US, scan-2, unweighted, directed
## TODO

test_that("Issue 18 is resolved", {

  library(igraph)
  g <- graph(c(1,2,2,1, 1,3,3,1, 2,4, 3,4, 3,5,5,3, 4,5,5,4))
  expect_that(local_scan(g, mode="all"), equals(c(4, 3, 7, 6, 5)))
  expect_that(local_scan(g, mode="out"), equals(c(4, 3, 7, 2, 5)))
  expect_that(local_scan(g, mode="in"), equals(c(4, 2, 4, 6, 5)))
})

test_that("Issue 18 is really resolved", {
  library(igraph)
  el <- c(1, 5, 1, 7, 2, 5, 2, 7, 2, 10, 2, 13, 2, 18, 3, 5, 3, 10, 3, 
          13, 4, 5, 4, 10, 5, 7, 5, 10, 5, 13, 5, 18, 6, 3, 6, 5, 6, 7, 
          6, 13, 7, 5, 8, 5, 8, 10, 8, 18, 9, 3, 9, 5, 9, 7, 9, 10, 11, 
          5, 12, 5, 12, 7, 14, 5, 14, 7, 14, 13, 14, 18, 15, 5, 15, 13, 
          15, 18, 16, 5, 16, 10, 16, 13, 16, 18, 17, 5)
  
  g <- graph(el)

  sc1 <- sapply(graph.neighborhood(g, order=1, mode="all"), ecount)
  sc2 <- local_scan(graph.us=g, mode="all", k=1)
  expect_that(sc1, equals(sc2))

  g2 <- induced.subgraph(g, 5:8)
  sc21 <- sapply(graph.neighborhood(g2, order=1, mode="all"), ecount)
  sc22 <- local_scan(graph.us=g2, mode="all", k=1)
  expect_that(sc21, equals(sc22))
})

test_that("Issue 20 is resolved", {

  library(igraph)
  set.seed(12345)
  g1 <- erdos.renyi.game(n=20, p=0.1, directed=TRUE)
  g2 <- erdos.renyi.game(n=20, p=0.1, directed=TRUE)
  ls <- local_scan(g2, g1, k=1, mode="all")
  correct <- c(4, 1, 2, 1, 1, 8, 1, 2, 0, 5, 2, 3, 3, 4, 5, 3, 5, 4, 2, 1)
  expect_that(ls, equals(correct))
})

test_that("FUN argument works, #32", {
  library(igraph)
  r1 <- local_scan(graph.ring(10), k=1, FUN="ecount")
  r2 <- local_scan(graph.ring(10), k=1, FUN=ecount)
  expect_that(r1, equals(rep(2, 10)))
  expect_that(r2, equals(rep(2, 10)))
})

