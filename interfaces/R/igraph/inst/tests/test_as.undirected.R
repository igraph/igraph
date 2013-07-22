
context("as.undirected")

test_that("as.undirected keeps attributes", {
  library(igraph)
  g <- graph.formula(A+-+B, A--+C, C+-+D)
  g$name <- "Tiny graph"
  E(g)$weight <- seq_len(ecount(g))

  g2 <- as.undirected(g, mode="collapse") ; df2 <- get.data.frame(g2)
  g3 <- as.undirected(g, mode="each")     ; df3 <- get.data.frame(g3)
  g4 <- as.undirected(g, mode="mutual")   ; df4 <- get.data.frame(g4)

  expect_that(g2$name, equals(g$name))
  expect_that(g3$name, equals(g$name))
  expect_that(g4$name, equals(g$name))
  
  expect_that(df2[order(df2[,1], df2[,2]),]$weight, equals(c(4,2,9)))
  expect_that(df3[order(df3[,1], df3[,2]),]$weight, equals(c(1,3,2,4,5)))
  expect_that(df4[order(df4[,1], df4[,2]),]$weight, equals(c(4,9)))
})
