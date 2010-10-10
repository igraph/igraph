
# R

library(igraph)

arpack2 <- function (func, extra = NULL, sym = FALSE,
                    options = igraph.arpack.default, 
                    env = parent.frame(), complex = !sym) {
  options.tmp <- igraph.arpack.default
  options.tmp[names(options)] <- options
  options <- options.tmp
  if (sym && complex) {
    complex <- FALSE
    warning("Symmetric matrix, setting `complex' to FALSE")
    }
  on.exit(.Call("R_igraph_finalizer", PACKAGE = "igraph"))
  .Call("R_igraph_arpack", func, extra, options, env, 
        sym, PACKAGE = "igraph")
}


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

arp2 <- function(M, sym=TRUE, which="LM", nev=1, maxiter=300, ...) {
  res <- arpack2(callback, extra=M, sym=sym,
                 options=list(n=nrow(M), nev=nev, ncv=2*nev+2,
                   which=which, maxiter=maxiter, ...))
  res$vectors <- cbind(res$vectors)
  res
}  

is.principal <- function(M, lambda, eps=1e-6) {
  e1 <- eigen(M)$values
  max(abs(e1)) - abs(lambda) < eps && min(abs(e1-lambda)) < eps
}
  
M <- matrix(sample(0:1, 100, replace=TRUE, prob=c(8,2)), 10)
res <- arp(M, nev=3, ncv=7, sym=FALSE, maxiter=10000)
res2 <- arp2(M, nev=3, ncv=7, sym=FALSE, maxiter=10000)

gc()

