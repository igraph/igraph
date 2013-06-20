
context("dyad.census")

test_that("dyad.census works", {
  library(igraph)
  ce <- simplify(read.graph(gzfile("celegansneural.gml.gz"), format="gml"))
  dc <- dyad.census(ce)

  expect_that(dc, equals(list(mut=197, asym=1951, null=41808)))
  expect_that(sum(is.mutual(ce)), equals(dc$mut * 2))
  expect_that(ecount(as.undirected(ce, mode="collapse")) - dc$mut,
              equals(dc$asym))
  expect_that(sum(unlist(dc)), equals(vcount(ce) * (vcount(ce)-1) / 2))
})
