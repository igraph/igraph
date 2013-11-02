
context("average.path.length")

test_that("average.path.length works", {
  library(igraph)

  apl <- function(graph) {
    sp <- shortest.paths(graph, mode="out")
    if (is.directed(graph)) {
      diag(sp) <- NA
    } else {
      sp[lower.tri(sp, diag=TRUE)] <- NA
    }
    sp[sp=="Inf"] <- NA
    mean(sp, na.rm=TRUE)
  }

  giant.component <- function(graph, mode="weak") {
    clu <- clusters(graph, mode=mode)
    induced.subgraph(graph, which(clu$membership==which.max(clu$csize)))
  }
  
  g <- giant.component(erdos.renyi.game(100, 3/100))
  expect_that(apl(g), equals(average.path.length(g)))

  g <- giant.component(erdos.renyi.game(100, 6/100, dir=TRUE), mode="strong")
  expect_that(apl(g), equals(average.path.length(g)))

  g <- erdos.renyi.game(100, 2/100)
  expect_that(apl(g), equals(average.path.length(g)))
  
  g <- erdos.renyi.game(100, 4/100, dir=TRUE)
  expect_that(apl(g), equals(average.path.length(g)))
})
