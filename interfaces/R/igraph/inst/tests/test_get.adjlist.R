
context("get.adjlist")

test_that("get.adjist works", {

  library(igraph)

  g <- erdos.renyi.game(50, 2/50)
  al <- get.adjlist(g)
  g2 <- graph.adjlist(al, mode="all")
  expect_that(graph.isomorphic(g, g2), is_true())
  expect_that(graph.isomorphic.vf2(g, g2, vertex.color1=1:vcount(g),
                                   vertex.color2=1:vcount(g2))$iso,
              is_true())

####

  el <- get.adjedgelist(g)
  for (i in 1:vcount(g)) {
    a <- as.numeric(E(g)[adj(i)])
    expect_that(length(a), equals(length(el[[i]])))
    expect_that(sort(el[[i]]), equals(sort(a)))
  }

  g <- erdos.renyi.game(50, 4/50, directed=TRUE)
  el1 <- get.adjedgelist(g, mode="out")
  el2 <- get.adjedgelist(g, mode="in")
  for (i in 1:vcount(g)) {
    a <- as.numeric(E(g)[from(i)])
    expect_that(length(a), equals(length(el1[[i]])))
    expect_that(sort(el1[[i]]), equals(sort(a)))
  }
  for (i in 1:vcount(g)) {
    a <- as.numeric(E(g)[to(i)])
    expect_that(length(a), equals(length(el2[[i]])))
    expect_that(sort(el2[[i]]), equals(sort(a)))
  }
  
})
