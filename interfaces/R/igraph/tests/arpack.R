
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

# TODO: further tests for typically hard cases
