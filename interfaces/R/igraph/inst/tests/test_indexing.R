
context("Indexing")

mm <- function(...) {
  v <- as.numeric(as.vector(list(...)))
  matrix(v, nrow=sqrt(length(v)))
}
am <- function(x) {
  x <- as.matrix(x)
  dimnames(x) <- NULL
  x
}

library(igraph)
library(Matrix, quietly=TRUE, warn.conflicts=FALSE)

g <- graph.tree(20)

test_that("[ indexing works", {
  
  ## Are these vertices connected?
  expect_that(g[1,2], equals(1))
  expect_that(am(g[c(1,1,7), c(2,3,14)]), equals(mm(1,1,0, 1,1,0, 0,0,1)))
  expect_that(am(g[c(1,1,7), c(5,3,12)]), equals(mm(0,0,0, 1,1,0 ,0,0,0)))
  expect_that(am(g[c(1,1,1,1), c(2,3,2,2)]), equals(matrix(1, 4, 4)))
  expect_that(am(g[c(8,17), c(17,8)]), equals(mm(1,0, 0,0)))

})

V(g)$name <- letters[1:vcount(g)]

test_that("[ indexing works with symbolic names", {

  ## The same with symbolic names
  expect_that(g['a','b'], equals(1))
  expect_that(am(g[c('a','a','g'), c('b','c','n')]),
              equals(mm(1,1,0, 1,1,0, 0,0,1)))
  expect_that(am(g[c('a','a','g'), c('e','c','l')]),
              equals(mm(0,0,0, 1,1,0, 0,0,0)))
  expect_that(am(g[c('a','a','a','a'), c('b','c','b','b')]),
              equals(matrix(1, 4, 4)))
  expect_that(am(g[c('h','q'), c('q','h')]), equals(mm(1,0, 0,0)))
})

test_that("[ indexing works with logical vectors", {

  ## Logical vectors
  lres <- structure(c(0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0), .Dim = c(2L, 20L),
                    .Dimnames = list(c("b", "c"), c("a", "b", "c",
                      "d", "e", "f", "g", "h", "i", "j", "k", "l",
                      "m", "n", "o", "p", "q", "r", "s", "t")))
  expect_that(g[degree(g,mode="in")==0,2], equals(1))
  expect_that(as.matrix(g[2:3,TRUE]), equals(lres))
})

test_that("[ indexing works with negative indices", {
  
  ## Negative indices
  nres <- structure(c(0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0), .Dim = c(2L, 19L),
                    .Dimnames=list(c("b", "c"),
                      c("b", "c", "d", "e", "f", "g", "h", "i", "j",
                        "k", "l", "m", "n", "o", "p", "q", "r", "s",
                        "t")))
  expect_that(as.matrix(g[2:3,-1]), equals(nres))
})

el <- get.edgelist(g, names=FALSE)
E(g)$weight <- el[,1] * el[,2]

test_that("[ indexing works with weighted graphs", {
  
  ## Weighted graphs
  expect_that(g[1,2], equals(2))
  expect_that(am(g[c(1,1,7), c(2,3,14)]), equals(mm(2,2,0, 3,3,0, 0,0,98)))
  expect_that(am(g[c(1,1,7), c(5,3,12)]), equals(mm(0,0,0, 3,3,0, 0,0,0)))
  expect_that(am(g[c(1,1,1,1), c(2,3,2,2)]),
              equals(mm(2,2,2,2, 3,3,3,3, 2,2,2,2, 2,2,2,2)))
  expect_that(am(g[c(8,17), c(17,8)]), equals(mm(136,0, 0,0)))
})

test_that("[ indexing works with weighted graphs and symbolic names",
          {
            
  ## Weighted graph, with symbolic names
  expect_that(g['a','b'], equals(2))
  expect_that(am(g[c('a','a','g'), c('b','c','n')]),
              equals(mm(2,2,0, 3,3,0, 0,0,98)))
  expect_that(am(g[c('a','a','g'), c('e','c','l')]),
              equals(mm(0,0,0, 3,3,0, 0,0,0)))
  expect_that(am(g[c('a','a','a','a'), c('b','c','b','b')]),
              equals(mm(2,2,2,2, 3,3,3,3, 2,2,2,2, 2,2,2,2)))  
  expect_that(am(g[c('h','q'), c('q','h')]), equals(mm(136,0, 0,0)))
})

################################################################

test_that("[[ indexing works", {

  ## Adjacent vertices
  expect_that(g[[1, ]], equals(list(a=2:3)))
  expect_that(g[[, 2]], equals(list(b=1)))
  expect_that(g[[, 2, directed=FALSE]], equals(list(b=c(1,4,5))))
  expect_that(g[[2, directed=FALSE]], equals(list(b=c(1,4,5))))

  expect_that(g[[1:3, ]], equals(list(a=2:3, b=4:5, c=6:7)))
  expect_that(g[[, 1:3]], equals(list(a=numeric(), b=1, c=1)))
})

test_that("[[ indexing works with symbolic names", {
  
  ## Same with vertex names
  expect_that(g[['a', ]], equals(list(a=2:3)))
  expect_that(g[[, 'b']], equals(list(b=1)))
  expect_that(g[[, 'b', directed=FALSE]], equals(list(b=c(1,4,5))))
  expect_that(g[['b', directed=FALSE]], equals(list(b=c(1,4,5))))

  expect_that(g[[letters[1:3],]], equals(list(a=2:3, b=4:5, c=6:7)))
  expect_that(g[[, letters[1:3]]], equals(list(a=numeric(), b=1, c=1)))
})

test_that("[[ indexing works with logical vectors", {

  ## Logical vectors
  expect_that(g[[degree(g,mode="in")==0,]], equals(list(a=2:3)))
})

test_that("[[ indexing works with filtering on both ends", {

  ## Filtering on both ends
  expect_that(g[[1:10, 1:10]], equals(list(a=2:3, b=4:5, c=6:7, d=8:9,
                                           e=10, f=numeric(),
                                           g=numeric(), h=numeric(),
                                           i=numeric(), j=numeric())))
})

################################################################

test_that("[ can query edge ids", {
  
  ## Query edge ids
  expect_that(g[1,2, edges=TRUE], equals(1))
  expect_that(am(g[c(1,1,7), c(2,3,14), edges=TRUE]),
              equals(mm(1,1,0, 2,2,0, 0,0,13)))
  expect_that(am(g[c(1,1,7), c(5,3,12), edges=TRUE]),
              equals(mm(0,0,0, 2,2,0, 0,0,0)))
  expect_that(am(g[c(1,1,1,1), c(2,3,2,2), edges=TRUE]),
              equals(mm(1,1,1,1, 2,2,2,2, 1,1,1,1, 1,1,1,1)))
  expect_that(am(g[c(8,17), c(17,8), edges=TRUE]),
              equals(mm(16,0, 0,0)))
})

test_that("[ can query edge ids with symbolic names", {
  
  ## The same with symbolic names
  expect_that(g['a','b', edges=TRUE], equals(1))
  expect_that(am(g[c('a','a','g'), c('b','c','n'), edges=TRUE]),
              equals(mm(1,1,0, 2,2,0, 0,0,13)))
  expect_that(am(g[c('a','a','g'), c('e','c','l'), edges=TRUE]),
              equals(mm(0,0,0, 2,2,0, 0,0,0)))
  expect_that(am(g[c('a','a','a','a'), c('b','c','b','b'), edges=TRUE]),
              equals(mm(1,1,1,1, 2,2,2,2, 1,1,1,1, 1,1,1,1)))
  expect_that(am(g[c('h','q'), c('q','h'), edges=TRUE]),
              equals(mm(16,0 ,0,0)))
})

################################################################

test_that("[[ can query incident edges", {
  
  ## Incident edges of vertices
  expect_that(g[[1, , edges=TRUE]], equals(list(a=1:2)))
  expect_that(g[[, 2, edges=TRUE]], equals(list(b=1)))
  expect_that(g[[, 2, directed=FALSE, edges=TRUE]],
              equals(list(b=c(3,4,1))))
  expect_that(g[[2, directed=FALSE, edges=TRUE]], equals(list(b=c(3,4,1))))

  expect_that(g[[1:3, , edges=TRUE]], equals(list(a=1:2, b=3:4, c=5:6)))
  expect_that(g[[, 1:3, edges=TRUE]], equals(list(a=numeric(), b=1, c=2)))
})

test_that("[[ queries edges with vertex names", {
  
  ## Same with vertex names
  expect_that(g[['a', , edges=TRUE]], equals(list(a=1:2)))
  expect_that(g[[, 'b', edges=TRUE]], equals(list(b=1)))
  expect_that(g[[, 'b', directed=FALSE, edges=TRUE]],
              equals(list(b=c(3,4,1))))
  expect_that(g[['b', directed=FALSE, edges=TRUE]],
              equals(list(b=c(3,4,1))))

  expect_that(g[[letters[1:3],, edges=TRUE]],
              equals(list(a=1:2, b=3:4, c=5:6)))
  expect_that(g[[, letters[1:3], edges=TRUE]],
              equals(list(a=numeric(), b=1, c=2)))

  ## Filtering on both ends
  expect_that(g[[1:10, 1:10, edges=TRUE]],
              equals(list(1:2, 3:4, 5:6, 7:8, 9, numeric(), numeric(),
                          numeric(), numeric(), numeric())))
})

#################################################################

test_that("[ handles from and to properly", {
  
  ## from & to
  g <- graph.tree(20)
  expect_that(g[from=c(1,2,2,3), to=c(3,4,8,7)], equals(c(1,1,0,1)))

  V(g)$name <- letters[1:20]
  expect_that(g[from=c("a","b","b","c"), to=c("c","d","h","g")],
              equals(c(1,1,0,1)))
  
  E(g)$weight <- (1:ecount(g)) ^ 2 
  expect_that(g[from=c("a","b","b","c"), to=c("c","d","h","g")],
              equals(c(4,9,0,36)))

  expect_that(g[from=c("a","b","b","c"), to=c("c","d","h","g"),
                edges=TRUE], equals(c(2,3,0,6)))
              

})
