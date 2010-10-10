
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

is.principal <- function(M, lambda, eps=1e-4) {
  e1 <- eigen(M)$values
  abs(max(abs(e1)) - abs(lambda)) < eps
}

for (i in 1:10000) {
  cat(".")
  g <- erdos.renyi.game(30, 1/30, directed=TRUE)
  M <- get.adjacency(g)
  clu <- clusters(g, mode="strong")$membership
  ae <- numeric()
  for (x in seq_along(unique(clu))) {
    v <- which(clu==x)
    if (length(v)==1) {
      ae <- c(ae, 0)
    } else if (length(v)<=13) {
      ae <- c(ae, eigen(M[v,v])$values)
    } else {
      M2 <- M[v,v]
      ae <- c(ae, arp(M2, nev=6, ncv=min(13,length(v)), sym=FALSE)$values)
    }
  }
  pae <- ae[ which.max(abs(ae)) ]
  if (!is.principal(M, pae)) { print("Gebasz2!"); break; }  
}
cat("\n")


