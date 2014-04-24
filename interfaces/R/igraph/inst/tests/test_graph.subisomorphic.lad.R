
context("graph.subisomorphic.lad")

test_that("graph.subisomorphic.lad works", {

  library(igraph)

  pattern <- graph.formula(1:2:3:4:5,
                           1 - 2:5, 2 - 1:5:3, 3 - 2:4, 4 - 3:5, 5 - 4:2:1)
  target <- graph.formula(1:2:3:4:5:6:7:8:9,
                          1 - 2:5:7, 2 - 1:5:3, 3 - 2:4, 4 - 3:5:6:8:9,
                          5 - 1:2:4:6:7, 6 - 7:5:4:9, 7 - 1:5:6,
                          8 - 4:9, 9 - 6:4:8)
  domains <- list(`1` = c(1,3,9), `2` = c(5,6,7,8), `3` = c(2,4,6,7,8,9),
                  `4` = c(1,3,9), `5` = c(2,4,8,9))
  i1 <- graph.subisomorphic.lad(pattern, target, all.maps=TRUE)
  i2 <- graph.subisomorphic.lad(pattern, target, induced=TRUE,
                                all.maps=TRUE)
  i3 <- graph.subisomorphic.lad(pattern, target, domains=domains,
                                all.maps=TRUE)

  expect_that(i1$iso, is_true())
  expect_that(i2, equals(
    structure(list(iso = TRUE, map = structure(c(1, 2, 3, 4, 5),
                                 .Names = c("1", "2", "3", "4", "5")),
                   maps = list(structure(c(1, 2, 3, 4, 5),
                     .Names = c("1", "2", "3", "4", "5")),
                     structure(c(6, 4, 3, 2, 5),
                               .Names = c("6", "4", "3", "2", "5")),
                     structure(c(6, 5, 2, 3, 4),
                               .Names = c("6", "5", "2", "3", "4")),
                     structure(c(1, 5, 4, 3, 2),
                               .Names = c("1", "5", "4", "3", "2")))),
              .Names = c("iso", "map", "maps")) ))
  expect_that(i3, equals(
    structure(list(iso = TRUE, map = structure(c(1, 5, 4, 3, 2),
                                 .Names = c("1", "5", "4", "3", "2")),
                   maps = list(structure(c(1, 5, 4, 3, 2),
                     .Names = c("1", "5", "4", "3", "2")))),
              .Names = c("iso", "map", "maps")) ))

})

test_that("LAD stress test", {

  library(igraph)
  set.seed(42)
  N <- 100
  
  for (i in 1:N) {
    target <- erdos.renyi.game(20, .5)
    pn <- sample(4:18, 1)
    pattern <- induced.subgraph(target, sample(vcount(target), pn))
    iso <- graph.subisomorphic.lad(pattern, target, induced=TRUE,
                                   all.maps=FALSE)
    expect_that(iso$iso, is_true())
  }

  set.seed(42)
  
  for (i in 1:N) {
    target <- erdos.renyi.game(20, 1/20)
    pn <- sample(5:18, 1)
    pattern <- erdos.renyi.game(pn, .6)
    iso <- graph.subisomorphic.lad(pattern, target, induced=TRUE,
                                   all.maps=FALSE)
    expect_that(iso$iso, is_false())
  }

})
