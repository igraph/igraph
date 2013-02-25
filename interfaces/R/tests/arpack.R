
library(igraph)
library(Matrix, warn.conflicts=FALSE)

####

f <- function(x, extra=NULL) x
res <- arpack(f, options=list(n=10, nev=2, ncv=4), sym=TRUE)
res$values

####

f <- function(x, extra=NULL) {
  y <- x
  y[1] <- (length(x)-1)*x[1] - sum(x[-1])
  for (i in 2:length(x)) {
    y[i] <- x[i] - x[1]
  }
  y
}
r1 <- arpack(f, options=list(n=10, nev=1, ncv=3), sym=TRUE)
r2 <- eigen(graph.laplacian(graph.star(10, mode="undirected")))

correctSign <- function(x) { if (x[1]<0) { -x } else { x } }
max(abs(r1$values - r2$values[1])) < 1e-14
max(abs(correctSign(r1$vectors) - correctSign(r2$vectors[,1]))) < 1e-14

####
# Complex case

A <- structure(c(-6, -6, 7, 6, 1, -9, -3, 2, -9, -7, 0, 1, -7, 8, 
-7, 10, 0, 0, 1, 1, 10, 0, 8, -4, -4, -5, 8, 9, -6, 9, 3, 8, 
6, -1, 9, -9, -6, -3, -1, -7, 8, -4, -4, 10, 0, 5, -2, 0, 7, 
10, 1, 4, -8, 3, 5, 3, -7, -9, 10, -1, -4, -7, -1, 7, 5, -5, 
1, -4, 9, -2, 10, 1, -7, 7, 6, 7, -3, 0, 9, -5, -8, 1, -3, -3, 
-8, -7, -8, 10, 8, 7, 0, 6, -7, -8, 10, 10, 1, 0, -2, 6), .Dim = c(10L, 
10L))

f <- function(x, extra=NULL) A %*% x
res <- arpack(f, options=list(n=10, nev=3, ncv=7, sym=FALSE))
res$values
(res$values[1] * res$vectors[,1]) / (A %*% res$vectors[,1])
(res$values[2] * res$vectors[,2]) / (A %*% res$vectors[,2])
abs((res$values[3] * res$vectors[,3]) / (A %*% res$vectors[,3]))

f <- function(x, extra=NULL) A %*% x
res <- arpack(f, options=list(n=10, nev=4, ncv=9, sym=FALSE))
res$values
(res$values[1] * res$vectors[,1]) / (A %*% res$vectors[,1])
(res$values[2] * res$vectors[,2]) / (A %*% res$vectors[,2])
abs((res$values[3] * res$vectors[,3]) / (A %*% res$vectors[,3]))
abs((res$values[4] * res$vectors[,4]) / (A %*% res$vectors[,4]))

####

# TODO: further tests for typically hard cases
