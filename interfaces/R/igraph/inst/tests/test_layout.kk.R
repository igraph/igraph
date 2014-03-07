
context("Kamada-Kawai layouts")

test_that("Kamada-Kawai layout generator works", {

  library(igraph)
  g <- graph.ring(10)
  l <- layout.kamada.kawai(g, maxiter=50)
  if (Sys.info()["sysname"] == "Darwin") {
    expect_that(sum(l), equals(-1.13071769106689))
  } else {
    expect_that(sum(l), equals(-6.77278645472984e-05))
  }

  g <- graph.star(30)
  l <- layout.kamada.kawai(g, maxiter=500)
  if (Sys.info()["sysname"] == "Darwin") {
    expect_that(sum(l), equals(-85.6883999492408))
  } else {
    expect_that(sum(l), equals(-86.1405864709501))
  }

  g <- graph.ring(10)
  E(g)$weight <- rep(1:2, length.out=ecount(g))
  l <- layout.kamada.kawai(g, maxiter=500)
  if (Sys.info()["sysname"] == "Darwin") {
    expect_that(sum(l), equals(1.61069099387368))
  } else {
    expect_that(sum(l), equals(-1.83036635516248))
  }
  
})

test_that("3D Kamada-Kawai layout generator works", {

  library(igraph)
  g <- graph.star(30)
  l <- layout.kamada.kawai(g, maxiter=5000, dim=3)
  expect_that(sum(l), equals(61.0559727551764))

})
