
context("edge names")

test_that("edge names work", {

  library(igraph) 

  ## named edges
  igraph.options(print.edge.attributes = TRUE)
  g <- graph.ring(10)
  E(g)$name <- letters[1:ecount(g)]
  g2 <- delete.edges(g, c("b", "d", "e"))
  expect_that(get.edgelist(g2),
              equals(structure(c(1, 3, 6, 7, 8, 9, 1, 2, 4, 7, 8, 9,
                                 10, 10), .Dim = c(7L, 2L))))

  ## named vertices
  g <- graph.ring(10)
  V(g)$name <- letters[1:vcount(g)]
  g3 <- delete.edges(g, c("a|b", "f|g", "c|b"))
  expect_that(get.edgelist(g3),
              equals(structure(c("c", "d", "e", "g", "h", "i", "a",
                                 "d", "e", "f", "h", "i", "j", "j"),
                               .Dim = c(7L, 2L))))


  ## no names at all, but select edges based on vertices
  g <- graph.ring(10)
  g4 <- delete.edges(g, c("1|2", "8|7", "1|10"))
  expect_that(get.edgelist(g4),
              equals(structure(c(2, 3, 4, 5, 6, 8, 9, 3, 4, 5, 6, 7,
                                 9, 10), .Dim = c(7L, 2L))))


  ## mix edge names and vertex names
  g <- graph.ring(10)
  V(g)$name <- letters[1:vcount(g)]
  E(g)$name <- LETTERS[1:ecount(g)]
  g5 <- delete.edges(g, c("a|b", "F", "j|i"))
  expect_that(get.edgelist(g5),
              equals(structure(c("b", "c", "d", "e", "g", "h", "a",
                                 "c", "d", "e", "f", "h", "i", "j"),
                               .Dim = c(7L, 2L))))
})
