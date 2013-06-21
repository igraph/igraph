
context("get.all.shortest.paths")

test_that("get.all.shortest.paths works", {

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

  sortlist <- function(list) {
    list <- lapply(list, sort)
    list[order(sapply(list, paste, collapse="!"))]
  }

  sp1 <- get.all.shortest.paths(g, "s", "t", weights=NA)

  expect_that(sortlist(sp1$res),
              equals(list(c(1, 2, 7), c(1, 3, 7))))
  expect_that(sp1$nrgeo,
              equals(c(1,1,1,1,1,1,2)))

  sp2 <- get.all.shortest.paths(g, "s", "t")

  expect_that(sortlist(sp2$res),
              equals(list(c(1, 2, 3, 4, 7), c(1, 2, 7), c(1, 3, 7))))
  expect_that(sp2$nrgeo, equals(c(1,1,2,1,1,1,3)))

  ## TODO

  ## E(g)$weight <- E(g)$weight - 1
  ## get.all.shortest.paths(g, "s", "t")

})
