
context("alpha.centrality")

test_that("dense alpha.centrality works", {
  library(igraph)
  g.1 <- graph( c(1,3,2,3,3,4,4,5) )
  ac1 <- alpha.centrality(g.1, sparse=FALSE)
  expect_that(ac1, equals(c(1, 1, 3, 4, 5)))
  
  g.2 <- graph( c(2,1,3,1,4,1,5,1) )
  ac2 <- alpha.centrality(g.2, sparse=FALSE)
  expect_that(ac2, equals(c(5,1,1,1,1)))
  
  g.3 <- graph( c(1,2,2,3,3,4,4,1,5,1) )
  ac3 <- alpha.centrality(g.3, alpha=0.5, sparse=FALSE)
  expect_that(ac3, equals(c(76, 68, 64, 62, 30)/30))
})

test_that("sparse alpha.centrality works", {
  if (require(Matrix, quietly=TRUE)) {
    library(igraph)
    g.1 <- graph( c(1,3,2,3,3,4,4,5) )
    ac1 <- alpha.centrality(g.1, sparse=TRUE)
    expect_that(ac1, equals(c(1, 1, 3, 4, 5)))
  
    g.2 <- graph( c(2,1,3,1,4,1,5,1) )
    ac2 <- alpha.centrality(g.2, sparse=TRUE)
    expect_that(ac2, equals(c(5,1,1,1,1)))
    
    g.3 <- graph( c(1,2,2,3,3,4,4,1,5,1) )
    ac3 <- alpha.centrality(g.3, alpha=0.5, sparse=TRUE)
    expect_that(ac3, equals(c(76, 68, 64, 62, 30)/30))
  }
})

##############################
## weighted version

test_that("weighted dense alpha.centrality works", {
  library(igraph)
  star <- graph.star(10)
  E(star)$weight <- sample(ecount(star))

  ac1 <- alpha.centrality(star, sparse=FALSE)
  expect_that(ac1, equals(c(46, 1, 1, 1, 1, 1, 1, 1, 1, 1)))

  ac2 <- alpha.centrality(star, weights="weight", sparse=FALSE)
  expect_that(ac2, equals(c(46, 1, 1, 1, 1, 1, 1, 1, 1, 1)))
  
  ac3 <- alpha.centrality(star, weights=NA, sparse=FALSE)
  expect_that(ac3, equals(c(vcount(star), 1, 1, 1, 1, 1, 1, 1, 1, 1)))
})

test_that("weighted sparse alpha.centrality works", {
  if (require("Matrix", quietly=TRUE)) {
    library(igraph)
    star <- graph.star(10)
    E(star)$weight <- sample(ecount(star))
    
    ac1 <- alpha.centrality(star, sparse=TRUE)
    expect_that(ac1, equals(c(46, 1, 1, 1, 1, 1, 1, 1, 1, 1)))
    
    ac2 <- alpha.centrality(star, weights="weight", sparse=TRUE)
    expect_that(ac2, equals(c(46, 1, 1, 1, 1, 1, 1, 1, 1, 1)))
    
    ac3 <- alpha.centrality(star, weights=NA, sparse=TRUE)
    expect_that(ac3, equals(c(vcount(star), 1, 1, 1, 1, 1, 1, 1, 1, 1)))
  }
})
