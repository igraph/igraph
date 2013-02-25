
library(igraph)

set.seed(42)                            # Shouldn't be here....

ashs <- function(graph, as=TRUE) {
  mscale <- function(x) {
    if (sd(x)!=0) { x <- scale(x) }
    if (x[1] < 0) { x <- -x       }
    x
  }
  A <- get.adjacency(graph, sparse=FALSE)
  if (as) { 
    s1 <- eigen(t(A) %*% A)$vectors[,1]
    s2 <- authority.score(graph)$vector
  } else {
    s1 <- eigen(A %*% t(A))$vectors[,1]
    s2 <- hub.score(graph)$vector
  }    
  all(abs(mscale(s1) - mscale(s2)) < 1e-12)
}

g1 <- ba.game(100, m=10)
ashs(g1)
ashs(g1, as=FALSE)

g2 <- erdos.renyi.game(100, 2/100)
ashs(g2)
ashs(g2, as=FALSE)

g3 <- graph.ring(100)
all(authority.score(g3)$vector == 1)
all(hub.score(g3)$vector == 1)

## Stress-test

is.principal <- function(M, lambda, eps=1e-12) {
  abs(eigen(M)$values[1] - lambda) < eps
}

is.ev <- function(M, v, lambda, eps=1e-12) {
  max(abs(M %*% v - lambda * v)) < eps
}

is.good <- function(M, v, lambda, eps=1e-12) {
  is.principal(M, lambda, eps) && is.ev(M, v, lambda, eps)
}

for (i in 1:1000) {
  G <- erdos.renyi.game(10, sample(1:20, 1), type="gnm")
  as <- authority.score(G)
  M <- get.adjacency(G)
  if (!is.good(t(M) %*% M, as$vector, as$value)) {
    print("Foobar!")
    break
  }
}

for (i in 1:1000) {
  G <- erdos.renyi.game(10, sample(1:20, 1), type="gnm")
  hs <- hub.score(G)
  M <- get.adjacency(G)
  if (!is.good(M %*% t(M), hs$vector, hs$value)) {
    print("Foobar!")
    break
  }
}
