
context("components")

test_that("components works", {
  library(igraph)
  set.seed(42)
  
  gc <- function(graph) {
    cl <- components(graph)
    induced_subgraph(graph, which(cl$membership==which.max(cl$csize)))
  }
  
  rg <- function(n) {
    gc(sample_gnp(n, 1/n))
  }
  
  G <- lapply(1:30, function(x) rg(sample(100, 1)))
  Gsize <- sapply(G, vcount)

  allg <- disjoint_union(G)
  clu <- components(allg)

  expect_that(as.numeric(table(clu$membership)), equals(clu$csize))
  expect_that(sort(clu$csize), equals(sort(Gsize)))
  expect_that(clu$no, equals(length(G)))
})

test_that("components names results", {
  library(igraph)

  g <- make_ring(10) + make_full_graph(5)
  V(g)$name <- letters[1:15]

  clu <- components(g)
  expect_that(names(clu$membership), equals(letters[1:15]))
})

test_that("groups works", {
  library(igraph)

  g <- make_ring(10) + make_full_graph(5)
  gr <- groups(components(g))

  expect_that(gr, equals(structure(list(`1` = 1:10, `2` = 11:15), .Dim = 2L,
                                   .Dimnames = list( c("1", "2")))))

  V(g)$name <- letters[1:15]
  gr <- groups(components(g))

  expect_that(gr, equals(structure(list(`1` = letters[1:10],
                                        `2` = letters[11:15]), .Dim = 2L,
                                        .Dimnames = list(c("1", "2")))))
})
