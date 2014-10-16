
context("Kamada-Kawai layouts")

test_that("Kamada-Kawai layout generator works", {

  library(igraph)
  g <- make_ring(10)
  l <- layout_with_kk(g, maxiter=50)
  if (Sys.info()["sysname"] == "Darwin") {
    expect_that(sum(l), equals(-1.13071769106689))
  } else if (Sys.info()["sysname"] == "Linux" &&
             Sys.info()["machine"] == "x86_64") {
    expect_that(sum(l), equals(-6.77278645472984e-05))
  } else if (Sys.info()["sysname"] == "Linux" &&
             Sys.info()["machine"] == "i686") {
    expect_that(sum(l), equals(0.914809637353466))
  }

  g <- make_star(30)
  l <- layout_with_kk(g, maxiter=500)
  if (Sys.info()["sysname"] == "Darwin") {
    expect_that(sum(l), equals(-85.6883999492408))
  } else if (Sys.info()["sysname"] == "Linux" &&
             Sys.info()["machine"] == "x86_64") {
    expect_that(sum(l), equals(-86.1405864709501))
  } else if (Sys.info()["sysname"] == "Linux" &&
             Sys.info()["machine"] == "i686") {
    expect_that(sum(l), equals(-85.142223229617))
  }

  g <- make_ring(10)
  E(g)$weight <- rep(1:2, length.out=ecount(g))
  l <- layout_with_kk(g, maxiter=500)
  if (Sys.info()["sysname"] == "Darwin") {
    expect_that(sum(l), equals(1.61069099387368))
  } else if (Sys.info()["sysname"] == "Linux" &&
             Sys.info()["machine"] == "x86_64") {
    expect_that(sum(l), equals(-1.83036635516248))
  } else if (Sys.info()["sysname"] == "Linux" &&
             Sys.info()["machine"] == "i686") {
    expect_that(sum(l), equals(0.0631144692360025))
  }
  
})

test_that("3D Kamada-Kawai layout generator works", {

  library(igraph)
  g <- make_star(30)
  l <- layout_with_kk(g, maxiter=5000, dim=3)
  expect_that(sum(l), equals(61.0559727551764))

})
