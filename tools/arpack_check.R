
# R

library(igraph)

callback <- function(x, extra) {
  y <- extra %*% x
}

residual <- function(A, v, lambda) {
  sum(abs( A %*% v - lambda * v )^2)
}

arp <- function(M, sym=TRUE, which="LM", nev=1, maxiter=300, ...) {
  res <- arpack(callback, extra=M, sym=sym,
                options=list(n=nrow(M), nev=nev, ncv=2*nev+2,
                  which=which, maxiter=maxiter, ...))
  res$vectors <- cbind(res$vectors)
  res
}  

is.principal <- function(M, lambda, eps=1e-6) {
  e1 <- eigen(M)$values
  max(abs(e1)) - abs(lambda) < eps && min(abs(e1-lambda)) < eps
}
  
## Symmetric solver, dense matrix

for (i in 1:10000) {
  M <- matrix(sample(1:10, 100, replace=TRUE), 10)
  M <- M+t(M)
  res <- arp(M)
  if (residual(M, res$vectors, res$values) > 1e-12) {
    print("Gebasz!" ); break;
  }
  if (!is.principal(M, res$values))     { print("Gebasz2!"); break; }
}

## Symmetric solver, sparse matrix

for (i in 1:10000) {
  M <- matrix(sample(0:1, 100, replace=TRUE, prob=c(8,2)), 10)
  M <- M+t(M)
  res <- arp(M, nev=1, ncv=4)
  w <- which.max(abs(res$values))
  if (residual(M, res$vectors[,w], res$values[w]) > 1e-12) {
    print("Gebasz!" ); break;
  }
  if (!is.principal(M, res$values[w])) { print("Gebasz2!"); break; }
}

## Symmetric solver, even sparser matrix

for (i in 1:10000) {
  M <- matrix(sample(0:1, 100, replace=TRUE, prob=c(14,1)), 10)
  M <- M+t(M)
  res <- arp(M, nev=1, ncv=5, maxiter=1e7)
  w <- which.max(abs(res$values))
  if (residual(M, res$vectors[,w], res$values[w]) > 1e-12) {
    print("Gebasz!" ); break;
  }
  if (!is.principal(M, res$values[w])) { print("Gebasz2!"); break; }
}

## Hard one:

M <- get.adjacency(graph.formula(a-b-c, d-e-f, g-h-i))
res <- arp(M, nev=2, ncv=5)
w <- which.max(abs(res$values))
residual(M, res$vectors[,w], res$values[w])
is.principal(M, res$values[w])

## Solution: nev=2, ncv=5, and always check which eigenvalue is
## larger.

#################################################################

## Non-symmetric case, dense

for (i in 1:10000) {
  M <- matrix(sample(1:10, 100, replace=TRUE), 10)
  res <- arp(M, sym=FALSE)
  if (residual(M, res$vectors, res$values) > 1e-12) {
    print("Gebasz!" ); break;
  }
  if (!is.principal(M, res$values))     { print("Gebasz2!"); break; }
}

## Non-symmetric, sparse

for (i in 1:10000) {
  M <- matrix(sample(0:1, 100, replace=TRUE, prob=c(8,2)), 10)
  res <- arp(M, nev=3, ncv=7, sym=FALSE)
  w <- which.max(abs(res$values))
  if (residual(M, res$vectors[,w], res$values[w]) > 1e-12) {
    print("Gebasz!" ); break;
  }
  if (!is.principal(M, res$values[w])) { print("Gebasz2!"); break; }
}

## Non-symmetric solver, even sparser matrix

for (i in 1:10000) {
  cat(".")
  M <- matrix(sample(0:1, 100, replace=TRUE, prob=c(9,1)), 10)
  res <- arp(M, nev=3, ncv=7, sym=FALSE, maxiter=10000)
  w <- which.max(abs(res$values))
  if (residual(M, res$vectors[,w], res$values[w]) > 1e-12) {
    print("Gebasz!" ); break;
  }
  if (!is.principal(M, res$values[w])) { print("Gebasz2!"); break; }
}

############################

M <- matrix(byrow=TRUE, nc=9,
            c(0,0,0,0,0,0,0,0,0,
              0,0,0,0,0,0,0,0,1,
              0,0,0,0,0,0,0,0,1,
              0,0,0,0,0,0,0,0,0,
              0,1,0,1,0,0,0,0,0,
              0,1,0,0,0,0,0,0,0,
              0,0,0,1,0,0,0,0,0,
              1,0,0,0,0,0,0,0,1,
              0,0,0,0,0,0,1,0,0))

arp(M, nev=3, ncv=7, sym=FALSE)$values


#########################
## Components

cc <- function(x) {
  g1 <- erdos.renyi.game(10, 2/10, directed=TRUE)
  g2 <- erdos.renyi.game(10, 2/10, directed=TRUE)
  g <- g1 %du% g2
  
  e1 <- eigen(get.adjacency(g1))$values
  e2 <- eigen(get.adjacency(g2))$values
  ee <- eigen(get.adjacency(g))$values

  max(abs(sort(Re(c(e1,e2))) - sort(Re(ee))),
      abs(sort(Im(c(e1,e2))) - sort(Im(ee))))
}

ccdiff <- sapply(1:1000, cc)


