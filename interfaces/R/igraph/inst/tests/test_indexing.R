
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

g <- make_tree(20)

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

el <- as_edgelist(g, names=FALSE)
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
  expect_that(g[[1, ]], is_equivalent_to(list(a=V(g)[2:3])))
  expect_that(g[[, 2]], is_equivalent_to(list(b=V(g)[1])))
  expect_that(g[[, 2, directed=FALSE]],
              is_equivalent_to(list(b=V(g)[c(1,4,5)])))
  expect_that(g[[2, directed=FALSE]],
              is_equivalent_to(list(b=V(g)[c(1,4,5)])))

  expect_that(g[[1:3, ]], is_equivalent_to(list(a=V(g)[2:3], b=V(g)[4:5],
                                                c=V(g)[6:7])))
  expect_that(g[[, 1:3]], is_equivalent_to(list(a=V(g)[numeric()],
                                                b=V(g)[1], c=V(g)[1])))
})

test_that("[[ indexing works with symbolic names", {
  
  ## Same with vertex names
  expect_that(g[['a', ]], is_equivalent_to(list(a=V(g)[2:3])))
  expect_that(g[[, 'b']], is_equivalent_to(list(b=V(g)[1])))
  expect_that(g[[, 'b', directed=FALSE]],
              is_equivalent_to(list(b=V(g)[c(1,4,5)])))
  expect_that(g[['b', directed=FALSE]],
              is_equivalent_to(list(b=V(g)[c(1,4,5)])))

  expect_that(g[[letters[1:3],]],
    is_equivalent_to(list(a=V(g)[2:3], b=V(g)[4:5], c=V(g)[6:7])))
  expect_that(g[[, letters[1:3]]],
    is_equivalent_to(list(a=V(g)[numeric()], b=V(g)[1], c=V(g)[1])))
})

test_that("[[ indexing works with logical vectors", {

  ## Logical vectors
  expect_that(g[[degree(g,mode="in")==0,]],
              is_equivalent_to(list(a=V(g)[2:3])))
})

test_that("[[ indexing works with filtering on both ends", {

  ## Filtering on both ends
  expect_that(g[[1:10, 1:10]],
    is_equivalent_to(list(a=V(g)[2:3], b=V(g)[4:5], c=V(g)[6:7], d=V(g)[8:9],
      e=V(g)[10], f=V(g)[numeric()], g=V(g)[numeric()], h=V(g)[numeric()],
                          i=V(g)[numeric()], j=V(g)[numeric()])))
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
  expect_that(g[[1, , edges=TRUE]], is_equivalent_to(list(a=E(g)[1:2])))
  expect_that(g[[, 2, edges=TRUE]], is_equivalent_to(list(b=E(g)[1])))
  expect_that(g[[, 2, directed=FALSE, edges=TRUE]],
              is_equivalent_to(list(b=E(g)[c(3,4,1)])))
  expect_that(g[[2, directed=FALSE, edges=TRUE]],
              is_equivalent_to(list(b=E(g)[c(3,4,1)])))

  expect_that(g[[1:3, , edges=TRUE]],
              is_equivalent_to(list(a=E(g)[1:2], b=E(g)[3:4], c=E(g)[5:6])))
  expect_that(g[[, 1:3, edges=TRUE]],
              is_equivalent_to(list(a=E(g)[numeric()], b=E(g)[1], c=E(g)[2])))
})

test_that("[[ queries edges with vertex names", {
  
  ## Same with vertex names
  expect_that(g[['a', , edges=TRUE]],
              is_equivalent_to(list(a=E(g)[1:2])))
  expect_that(g[[, 'b', edges=TRUE]],
              is_equivalent_to(list(b=E(g)[1])))
  expect_that(g[[, 'b', directed=FALSE, edges=TRUE]],
              is_equivalent_to(list(b=E(g)[c(3,4,1)])))
  expect_that(g[['b', directed=FALSE, edges=TRUE]],
              is_equivalent_to(list(b=E(g)[c(3,4,1)])))

  expect_that(g[[letters[1:3],, edges=TRUE]],
              is_equivalent_to(list(a=E(g)[1:2], b=E(g)[3:4], c=E(g)[5:6])))
  expect_that(g[[, letters[1:3], edges=TRUE]],
              is_equivalent_to(list(a=E(g)[numeric()], b=E(g)[1], c=E(g)[2])))

  ## Filtering on both ends
  expect_that(g[[1:10, 1:10, edges=TRUE]],
    is_equivalent_to(list(E(g)[1:2], E(g)[3:4], E(g)[5:6], E(g)[7:8],
                          E(g)[9], E(g)[numeric()], E(g)[numeric()],
                          E(g)[numeric()], E(g)[numeric()], E(g)[numeric()])))
})

#################################################################

test_that("[ handles from and to properly", {
  
  ## from & to
  g <- make_tree(20)
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

test_that("[[ works with from and to", {

  g <- make_tree(20)

  expect_equivalent(g[[1, ]], g[[from = 1]])
  expect_equivalent(g[[, 1]], g[[to = 1]])
  expect_equivalent(g[[1:5, 4:10]], g[[from = 1:5, to = 4:10]])

  expect_error(g[[1, from = 1]], "Cannot give both")
  expect_error(g[[, 2, to = 10]], "Cannot give both")
})

test_that("[[ returns vertex and edges sequences", {

  g <- make_tree(20)
  expect_true(is_igraph_vs(g[[1]][[1]]))
  expect_true(is_igraph_es(g[[1, edges = TRUE]][[1]]))
  expect_true(is_igraph_vs(g[[1:3, 2:6]][[1]]))
  expect_true(is_igraph_es(g[[1:3, 2:6, edges = TRUE]][[1]]))

})
