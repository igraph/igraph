
context("graph.subisomorphic, lad")

test_that("graph.subisomorphic, method = 'lad' works", {

  library(igraph)

  pattern <- graph_from_literal(1:2:3:4:5,
                       1 - 2:5, 2 - 1:5:3, 3 - 2:4, 4 - 3:5, 5 - 4:2:1)
  target <- graph_from_literal(1:2:3:4:5:6:7:8:9,
                      1 - 2:5:7, 2 - 1:5:3, 3 - 2:4, 4 - 3:5:6:8:9,
                      5 - 1:2:4:6:7, 6 - 7:5:4:9, 7 - 1:5:6,
                      8 - 4:9, 9 - 6:4:8)
  domains <- list(`1` = c(1,3,9), `2` = c(5,6,7,8), `3` = c(2,4,6,7,8,9),
                  `4` = c(1,3,9), `5` = c(2,4,8,9))
  i1 <- subgraph_isomorphic(pattern, target, method = "lad")
  i2 <- subgraph_isomorphic(pattern, target, induced=TRUE, method = "lad")
  i3 <- subgraph_isomorphic(pattern, target, domains=domains,
                            method = "lad")

  expect_that(i1, is_true())
  expect_that(i2, is_true())
  expect_that(i3, is_true())

})

test_that("LAD stress test", {

  library(igraph)
  set.seed(42)
  N <- 100
  
  for (i in 1:N) {
    target <- sample_gnp(20, .5)
    pn <- sample(4:18, 1)
    pattern <- induced_subgraph(target, sample(vcount(target), pn))
    iso <- subgraph_isomorphic(pattern, target, induced=TRUE,
                               method = "lad")
    expect_that(iso, is_true())
  }

  set.seed(42)
  
  for (i in 1:N) {
    target <- sample_gnp(20, 1/20)
    pn <- sample(5:18, 1)
    pattern <- sample_gnp(pn, .6)
    iso <- subgraph_isomorphic(pattern, target, induced=TRUE,
                               method = "lad")
    expect_that(iso, is_false())
  }

})
