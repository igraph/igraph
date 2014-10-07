
context("graph.adjancency")

test_that("graph_from_adjacency_matrix works", {

  library(igraph)

  M1 <- rbind(c(0,0,1,1),
              c(1,0,0,0),
              c(0,1,0,1),
              c(1,0,0,1))
  g1 <- graph_from_adjacency_matrix(M1)
  el1 <- as_edgelist(g1)
  expect_that(el1[order(el1[,1], el1[,2]),],
              equals(structure(c(1, 1, 2, 3, 3, 4, 4, 3, 4, 1, 2, 4,
                                 1, 4), .Dim = c(7L, 2L))))

  M2 <- rbind(c(0,1,1,1),
              c(1,0,0,0),
              c(1,0,0,1),
              c(1,0,1,0))
  g2 <- graph_from_adjacency_matrix(M2, mode="undirected")
  el2 <- as_edgelist(g2)
  expect_that(el2[order(el2[,1], el2[,2]),],
              equals(structure(c(1, 1, 1, 3, 2, 3, 4, 4),
                               .Dim = c(4L, 2L))))

  M3 <- rbind(c(0,1,1,2),
              c(1,0,0,0),
              c(1,0,0,0),
              c(1,0,1,0))
  g3 <- graph_from_adjacency_matrix(M3, mode="min")
  el3 <- as_edgelist(g3)
  expect_that(el3[order(el3[,1], el3[,2]),],
              equals(structure(c(1, 1, 1, 2, 3, 4), .Dim=c(3L, 2L))))

  M4 <- rbind(c(0,1,1,2),
              c(1,0,0,0),
              c(1,0,0,0),
              c(1,0,1,0))
  g4 <- graph_from_adjacency_matrix(M4, mode="max")
  el4 <- as_edgelist(g4)
  expect_that(el4[order(el4[,1], el4[,2]),],
              equals(structure(c(1, 1, 1, 1, 3,
                                 2, 3, 4, 4, 4), .Dim=c(5L, 2L))))

  M5 <- rbind(c(0,1,1,2),
              c(1,0,0,0),
              c(1,0,0,0),
              c(1,0,1,0))
  g5 <- graph_from_adjacency_matrix(M5, mode="upper")
  el5 <- as_edgelist(g5)
  expect_that(el5[order(el5[,1], el5[,2]),],
              equals(structure(c(1, 1, 1, 1, 2, 3, 4, 4), .Dim=c(4L, 2L))))

  M6 <- rbind(c(0,1,1,2),
              c(1,0,0,0),
              c(1,0,0,0),
              c(1,0,1,0))
  g6 <- graph_from_adjacency_matrix(M6, mode="lower")
  el6 <- as_edgelist(g6)
  expect_that(el6[order(el6[,1], el6[,2]),],
              equals(structure(c(1, 1, 1, 3, 2, 3, 4, 4), .Dim=c(4L, 2L))))

  M7 <- rbind(c(0,1,1,2),
              c(1,0,0,0),
              c(1,0,0,0),
              c(1,0,1,0))
  g7 <- graph_from_adjacency_matrix(M7, mode="plus")
  el7 <- as_edgelist(g7)
  expect_that(el7[order(el7[,1], el7[,2]),],
              equals(structure(c(1, 1, 1, 1, 1, 1, 1, 3, 2, 2, 3, 3,
                                 4, 4, 4, 4), .Dim = c(8L, 2L))))

  M8 <- rbind(c(0,1,1,0.5),
              c(1,0,0,0),
              c(1,0,0,0),
              c(1,0,2,0))
  g8 <- graph_from_adjacency_matrix(M8, mode="directed", weighted=TRUE)
  el8 <- cbind(as_edgelist(g8), E(g8)$weight)
  expect_that(el8[order(el8[,1], el8[,2]),],
              equals(structure(c(1, 1, 1, 2, 3, 4, 4, 2, 3, 4, 1, 1,
                                 1, 3, 1, 1, 0.5, 1, 1, 1, 2), .Dim =
                               c(7L, 3L))))

  M9 <- rbind(c(0,1,1,3),
              c(1,0,0,0),
              c(1,0,0,2),
              c(3,0,2,0))
  g9 <- graph_from_adjacency_matrix(M9, mode="undirected", weighted=TRUE)
  el9 <- cbind(as_edgelist(g9), E(g9)$weight)
  expect_that(el9[order(el9[,1], el9[,2]),],
              equals(structure(c(1, 1, 1, 3, 2, 3, 4, 4, 1, 1, 3, 2),
                               .Dim = c(4L, 3L))))

  M10 <- rbind(c(0,1,1,0.5),
               c(1,0,0,0),
               c(1,0,0,0),
               c(1,0,2,0))
  g10 <- graph_from_adjacency_matrix(M10, mode="max", weighted=TRUE)
  el10 <- cbind(as_edgelist(g10), E(g10)$weight)
  expect_that(el10[order(el10[,1], el10[,2]),],
              equals(structure(c(1, 1, 1, 3, 2, 3, 4, 4, 1, 1, 1, 2),
                               .Dim = c(4L, 3L))))

  M11 <- rbind(c(0,1,1,0.5),
               c(1,0,0,0),
               c(1,0,0,0),
               c(1,0,2,0))
  g11 <- graph_from_adjacency_matrix(M11, mode="min", weighted=TRUE)
  el11 <- cbind(as_edgelist(g11), E(g11)$weight)
  expect_that(el11[order(el11[,1], el11[,2]),],
              equals(structure(c(1, 1, 1, 2, 3, 4, 1, 1, 0.5), .Dim =
                               c(3L, 3L))))

  M12 <- rbind(c(0,1,1,0.5),
               c(1,0,0,0),
               c(1,0,0,0),
               c(1,0,2,0))
  g12 <- graph_from_adjacency_matrix(M12, mode="lower", weighted=TRUE)
  el12 <- cbind(as_edgelist(g12), E(g12)$weight)
  expect_that(el12[order(el12[,1], el12[,2]),],
              equals(structure(c(1, 1, 1, 3, 2, 3, 4, 4, 1, 1, 1, 2),
                               .Dim = c(4L, 3L))))

  M13 <- rbind(c(0,1,1,0.5),
               c(1,0,0,0),
               c(1,0,0,0),
               c(1,0,2,0))
  g13 <- graph_from_adjacency_matrix(M13, mode="upper", weighted=TRUE)
  el13 <- cbind(as_edgelist(g13), E(g13)$weight)
  expect_that(el13[order(el13[,1], el13[,2]),],
              equals(structure(c(1, 1, 1, 2, 3, 4, 1, 1, 0.5), .Dim =
                               c(3L, 3L))))

  M14 <- rbind(c(0,1,1,0.5),
               c(1,0,0,0),
               c(1,0,0,0),
               c(1,0,2,0))
  g14 <- graph_from_adjacency_matrix(M14, mode="plus", weighted=TRUE)
  el14 <- cbind(as_edgelist(g14), E(g14)$weight)
  expect_that(el14[order(el14[,1], el14[,2]),],
              equals(structure(c(1, 1, 1, 3, 2, 3, 4, 4, 2, 2, 1.5,
                                 2), .Dim = c(4L, 3L))))

})

test_that("graph_from_adjacency_matrix 2 edge bug is fixed", {

  library(Matrix)
  library(igraph)
  A <- Matrix(0, 10, 10, sparse=TRUE)
  A[3,5] <- A[5,3] <- 1
  g <- graph_from_adjacency_matrix(A, mode="undirected")
  expect_that(g[], equals(A))

})

test_that("graph.adjacenct empty graph bug is fixed", {

  library(Matrix)
  library(igraph)
  A <- Matrix(0, 10, 10, sparse=TRUE)
  g <- graph_from_adjacency_matrix(A, mode="undirected")
  expect_that(as.matrix(g[]), equals(as.matrix(A)))

})

test_that("bug #554 is fixed", {

  library(igraph)
  library(Matrix)

  M <- Matrix(0, 5, 5)
  M[1,2] <- M[2,1] <- M[3,4] <- M[4,3] <- 1
  g <- graph_from_adjacency_matrix(M, mode="undirected", weighted=TRUE)
  expect_that(g[], equals(M))

})
