
context("mean_distance")

test_that("mean_distance works", {
  library(igraph)

  apl <- function(graph) {
    sp <- distances(graph, mode="out")
    if (is_directed(graph)) {
      diag(sp) <- NA
    } else {
      sp[lower.tri(sp, diag=TRUE)] <- NA
    }
    sp[sp=="Inf"] <- NA
    mean(sp, na.rm=TRUE)
  }

  giant.component <- function(graph, mode="weak") {
    clu <- components(graph, mode=mode)
    induced_subgraph(graph, which(clu$membership==which.max(clu$csize)))
  }
  
  g <- giant.component(sample_gnp(100, 3/100))
  expect_that(apl(g), equals(mean_distance(g)))

  g <- giant.component(sample_gnp(100, 6/100, dir=TRUE), mode="strong")
  expect_that(apl(g), equals(mean_distance(g)))

  g <- sample_gnp(100, 2/100)
  expect_that(apl(g), equals(mean_distance(g)))
  
  g <- sample_gnp(100, 4/100, dir=TRUE)
  expect_that(apl(g), equals(mean_distance(g)))
})
