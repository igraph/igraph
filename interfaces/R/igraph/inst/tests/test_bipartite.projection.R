
context("bipartite.projection")

test_that("bipartite.projection works", {
  library(igraph)
  set.seed(42)

  g <- graph.full.bipartite(10,5)
  proj <- bipartite.projection(g)
  expect_that(graph.isomorphic(proj[[1]], graph.full(10)), is_true())
  expect_that(graph.isomorphic(proj[[2]], graph.full(5)), is_true())

  M <- matrix(0, nr=5, nc=3)
  rownames(M) <- c("Alice", "Bob", "Cecil", "Dan", "Ethel")
  colnames(M) <- c("Party", "Skiing", "Badminton")
  M[] <- sample(0:1, length(M), replace=TRUE)
  M
  g2 <- graph.incidence(M)
  expect_that(as.matrix(g2[1:5,6:8]), equals(M))
  expect_that(as.matrix(g2[1:5,1:5]), is_equivalent_to(matrix(0, 5, 5)))
  expect_that(as.matrix(g2[6:8,6:8]), is_equivalent_to(matrix(0, 3, 3)))  
    
  g2$name <- "Event network"
  proj2 <- bipartite.projection(g2)
  expect_that(as.matrix(proj2[[1]][]),
              is_equivalent_to(cbind(c(0,2,0,2,2), c(2,0,1,2,2),
                                     c(0,1,0,0,0), c(2,2,0,0,2),
                                     c(2,2,0,2,0))))
  expect_that(as.matrix(proj2[[2]][]),
              is_equivalent_to(cbind(c(0,4,1), c(4,0,1), c(1,1,0))))
  
  bs <- bipartite.projection.size(g2)
  expect_that(bs$vcount1, equals(vcount(proj2[[1]])))
  expect_that(bs$ecount1, equals(ecount(proj2[[1]])))
  expect_that(bs$vcount2, equals(vcount(proj2[[2]])))
  expect_that(bs$ecount2, equals(ecount(proj2[[2]])))
})

test_that("bipartite.projection can calculate only one projection", {
  library(igraph)
  set.seed(42)

  g <- bipartite.random.game(5, 10, p=.3)
  proj <- bipartite.projection(g)
  proj1 <- bipartite.projection(g, which="false")
  proj2 <- bipartite.projection(g, which="true")

  expect_that(graph.isomorphic(proj$proj1, proj1), is_true())
  expect_that(graph.isomorphic(proj$proj2, proj2), is_true())
  expect_that(vertex.attributes(proj$proj1), equals(vertex.attributes(proj1)))
  expect_that(vertex.attributes(proj$proj2), equals(vertex.attributes(proj2)))
  expect_that(edge.attributes(proj$proj1), equals(edge.attributes(proj1)))
  expect_that(edge.attributes(proj$proj2), equals(edge.attributes(proj2)))

})

