
context("GEM layout")

test_that("GEM layout works", {

  set.seed(42)
  l <- cbind(1:10, rep(2:3, each=5))
  g <- make_ring(10)
  l2 <- layout_with_gem(g, coords=l, maxiter=1)
  test_that(l[1:9,], equals(l2[1:9,]))  
  
})
