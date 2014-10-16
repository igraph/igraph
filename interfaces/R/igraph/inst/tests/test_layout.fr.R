
context("Fruchterman-Reingold layout")

test_that("", {

  library(igraph)
  set.seed(42)
  g <- make_ring(10)
  l <- layout_with_fr(g, niter=50, start.temp=sqrt(10)/10)
  if (.Machine$sizeof.pointer == 4) {
    expect_that(sum(l), equals(10.794223604849))
  } else {
    expect_that(sum(l), equals(10.7943032688805))
  }

  set.seed(42)
  g <- make_star(30)
  l <- layout_with_fr(g, niter=500, dim=3, start.temp=20)
  if (.Machine$sizeof.pointer == 4) {
    expect_that(sum(l), equals(1004.00737470853))
  } else {
    expect_that(sum(l), equals(941.472420651506))
  }

})
