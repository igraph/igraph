
context("leading.eigenvector.community")

test_that("leading.eigenvector.community works", {

  library(igraph)

  ## Check-test

  f <- function(membership, community, value, vector, multiplier, extra) {
    M <- sapply(1:length(vector), function(x) {
      v <- rep(0, length(vector))
      v[x] <- 1
      multiplier(v)
    })
    ev <- eigen(M)
    ret <- 0
    expect_that(ev$values[1], equals(value))
    if (sign(ev$vectors[1,1]) != sign(vector[1])) {
      ev$vectors <- -ev$vectors
    }
    expect_that(ev$vectors[,1], equals(vector))
  }

  g <- graph.famous("Zachary")
  lc <- leading.eigenvector.community(g, callback=f)
  
  expect_that(lc$modularity, equals(modularity(g, lc$membership)))
  expect_that(membership(lc),
              equals(c(1, 3, 3, 3, 1, 1, 1, 3, 2, 2, 1, 1, 3, 3, 2, 2,
                       1, 3, 2, 3, 2, 3, 2, 4, 4, 4, 2, 4, 4, 2, 2, 4,
                       2, 2)))
  expect_that(length(lc), equals(4))
  expect_that(sizes(lc),
              equals(structure(c(7L, 12L, 9L, 6L), .Dim = 4L, .Dimnames =
                               structure(list(`Community sizes` =
                                              c("1", "2", "3", "4")),
                                         .Names = "Community sizes"),
                               class = "table")))

  ## Check that the modularity matrix is correct

  f <- function(membership, community, value, vector, multiplier, extra) {
    M <- sapply(1:length(vector), function(x) {
      v <- rep(0, length(vector))
      v[x] <- 1
      multiplier(v)
    })
    myc <- membership==community
    B <- A[myc,myc] - (deg[myc] %*% t(deg[myc]))/2/ec
    BG <- B-diag(rowSums(B))
    
    expect_that(M, equals(BG))
  }

  g <- graph.famous("Zachary")
  A <- get.adjacency(g, sparse=FALSE)
  ec <- ecount(g)
  deg <- degree(g)
  lc <- leading.eigenvector.community(g, callback=f)

  ## Stress-test

  for (i in 1:100) {
    g <- erdos.renyi.game(20, sample(5:40, 1), type="gnm")
    lec1 <- leading.eigenvector.community(g)
    lec2 <- leading.eigenvector.community(g)
    expect_that(membership(lec1), equals(membership(lec2)))
  }

})
