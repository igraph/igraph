
context("operators")

test_that("operators work", {

  library(igraph)

  o <- function(x) x[order(x[,1], x[,2]),]

  g1 <- ring(10)
  g2 <- star(11, center=11, mode="undirected")
  gu <- union(g1, g2)
  expect_that(vcount(gu), equals(11))
  expect_that(ecount(gu), equals(20))
  expect_that(o(rbind(as_edgelist(g1), as_edgelist(g2))),
              equals(o(as_edgelist(gu))))

  gdu <- disjoint_union(g1, g2)
  expect_that(o(as_edgelist(gdu)),
              equals(o(rbind(as_edgelist(g1),
                             as_edgelist(g2)+vcount(g1)))))

####

  expect_that(graph.isomorphic(difference(gu, g1), g2), is_true())

####

  expect_that(graph.isomorphic(graph.intersection(gu, g2), g2), is_true())

  expect_that(graph.isomorphic(graph.intersection(gu, g1,
                                                  keep.all.vertices=FALSE),
                               g1),is_true())

####

  expect_that(complementer(complementer(g2)),
              equals(g2))

####

  gc <- compose(gu, g1)
  expect_that(vcount(gc), equals(11))
  expect_that(ecount(gc), equals(60))
  expect_that(diameter(gc), equals(2))
  
})
