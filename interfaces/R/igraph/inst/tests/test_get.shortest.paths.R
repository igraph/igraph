
context("get.shortest.paths")

test_that("get.shortest.paths works", {

  library(igraph)

  edges <- matrix(c("s", "a", 2,
                    "s", "b", 4,
                    "a", "t", 4,
                    "b", "t", 2,
                    "a", "1", 1,
                    "a", "2", 1,
                    "a", "3", 2,
                    "1", "b", 1,
                    "2", "b", 2,
                    "3", "b", 1),
                  byrow=TRUE, ncol=3,
                  dimnames=list(NULL, c("from", "to", "weight")))
  edges <- as.data.frame(edges)
  edges[[3]] <- as.numeric(as.character(edges[[3]]))

  g <- graph.data.frame(as.data.frame(edges))

  all1 <- get.all.shortest.paths(g, "s", "t", weights=NA)$res
  all2 <- get.all.shortest.paths(g, "s", "t")$res

  s1 <- get.shortest.paths(g, "s", "t", weights=NA)
  s2 <- get.shortest.paths(g, "s", "t")

  expect_that(s1$vpath %in% all1, is_true()) 
  expect_that(s2$vpath %in% all2, is_true())

})
