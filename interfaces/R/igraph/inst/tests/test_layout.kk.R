
context("Kamada-Kawai layouts")

test_that("Kamada-Kawai layout generator works", {

  library(igraph)
  g <- graph.ring(10)
  l <- layout.kamada.kawai(g, maxiter=50)
  expect_that(sum(l), equals(-1.13071769106689))

  g <- graph.star(30)
  l <- layout.kamada.kawai(g, maxiter=500)
  expect_that(sum(l), equals(-85.6883999492408))

  g <- graph.ring(10)
  E(g)$weight <- rep(1:2, length.out=ecount(g))
  l <- layout.kamada.kawai(g, maxiter=500)
  expect_that(sum(l), equals(1.61069099387368))
  
})

test_that("3D Kamada-Kawai layout generator works", {

  library(igraph)
  g <- graph.star(30)
  l <- layout.kamada.kawai(g, maxiter=5000, dim=3)
  expect_that(sum(l), equals(61.0559727551764))

})
