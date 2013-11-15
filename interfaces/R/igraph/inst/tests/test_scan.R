
context("Local scan statistics")

test_that("Local scan-1 statistics works", {
  
  library(igraph)
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

my.local.scan <- function(g, gp=NULL, k=1, mode="out", FUN=ecount,
                          weighted=FALSE) {

  wstat <- function(g,...) {
    A <- get.adjacency(g, attr="weight")
    if (k==0) { # weighted degree
      indeg <- colSums(A)
      outdeg <- rowSums(A)
      if (mode=="in") {
        out <- indeg
      } else if (mode=="out") {
        out <- outdeg
      } else {
        out <- indeg + outdeg
      }
    } else { # weighted ecount
      out <- sum(A)
      if (mode=="in" | mode=="out") out <- out/2
    }
    return(out)
  }

  if(k < 0) { stop("Error: k should be a non-negative integer!\n") }
  if(k == 0) { FUN <- degree }
  if(weighted) { FUN <- wstat }

  require(igraph)

  if (is.matrix(g) | is.matrix(gp)) {
    gmode <- ifelse((mode=="out" | mode=="in"),"directed","undirected")
    g <- simplify(graph.adjacency(g, mode=gmode))
    if (!is.null(gp)) { gp <- simplify(graph.adjacency(gp, mode=gmode)) }
  }

  n <- vcount(g)
  if (is.null(gp)) {
    if (k==0) {
      out <- FUN(g,mode=mode)
    } else {
      out <- sapply(graph.neighborhood(g, k, V(g), mode), FUN)
    }
  } else { # them
    if (k==0) {
      out <- unlist(sapply(1:n, function(x) {
        vid <- unlist(neighborhood(g, k+1, x, mode));
        x.rank <- which(sort(vid)==x);
        FUN(induced.subgraph(gp, vid), mode=mode)[x.rank]}))
    } else {
      out <- sapply(V(g), function(x) {
        FUN(induced.subgraph(gp,unlist(neighborhood(g, k, x, mode))))})
    }
  }
  return(out)
}

test_that("General scan-stat works", {

  require(igraph)
  require(Matrix)

  set.seed(12345)
  n <- 10^3
  p <- 0.1
  g <- erdos.renyi.game(n,p)
  E(g)$weight = sample(ecount(g))
  gp <- erdos.renyi.game(n,p)
  E(gp)$weight = sample(ecount(gp))

  ## us: scan0, unweighted
  system.time(s1 <- local.scan(g, k=0))
  system.time(s2 <- my.local.scan(g, k=0))
  expect_that(s1, equals(s2))

  ## us: scan0, weighted
  system.time(s1 <- local.scan(g, k=0, weighted=TRUE))
  system.time(s2 <- my.local.scan(g, k=0, weighted=TRUE))
  expect_that(s1, equals(s2))

  ## us: scan1, unweighted
  system.time(s1 <- local.scan(g))
  system.time(s2 <- my.local.scan(g, k=1))
  expect_that(s1, equals(s2))

  ## us: scan1, weighted
  system.time(s1 <- local.scan(g, k=1, weighted=TRUE))
  system.time(s2 <- my.local.scan(g, k=1, weighted=TRUE))
  expect_that(s1, equals(s2))

  ## us: scan2, unweighted
  system.time(s1 <- local.scan(g, k=2))
  system.time(s2 <- my.local.scan(g, k=2))
  expect_that(s1, equals(s2))

  ## us: scan2, weighted
  system.time(s1 <- local.scan(g, k=2, weighted=TRUE))
  system.time(s2 <- my.local.scan(g, k=2, weighted=TRUE))
  expect_that(s1, equals(s2))

  ## them: scan0, unweighted
  system.time(s1 <- local.scan(g, gp, k=0))
  system.time(s2 <- my.local.scan(g, gp, k=0))
  expect_that(s1, equals(s2))

  ## them: scan0, weighted
  system.time(s1 <- local.scan(g, gp, k=0, weighted=TRUE))
  system.time(s2 <- my.local.scan(g, gp, k=0, weighted=TRUE))
  expect_that(s1, equals(s2))

  ## them: local.scan, unweighted
  system.time(s1 <- local.scan(g, gp, k=1))
  system.time(s2 <- my.local.scan(g, gp, k=1))
  expect_that(s1, equals(s2))

  ## them: local.scan, weighted
  system.time(s1 <- local.scan(g, gp, k=1, weighted=TRUE))
  system.time(s2 <- my.local.scan(g, gp, k=1, weighted=TRUE))
  expect_that(s1, equals(s2))

  ## them: scan2, unweighted
  system.time(s1 <- local.scan(g, gp, k=2))
  system.time(s2 <- my.local.scan(g, gp, k=2))
  expect_that(s1, equals(s2))

  ## them: scan2, weighted
  system.time(s1 <- local.scan(g, gp, k=2, weighted=TRUE))
  system.time(s2 <- my.local.scan(g, gp, k=2, weighted=TRUE))
  expect_that(s1, equals(s2))

})
