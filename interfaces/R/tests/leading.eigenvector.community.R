
library(igraph)

## Check-test

f <- function(membership, community, value, vector, multiplier, extra) {
  M <- sapply(1:length(vector), function(x) {
    v <- rep(0, length(vector))
    v[x] <- 1
    multiplier(v)
  })
  ev <- eigen(M)
  ret <- 0
  if (abs(ev$values[1] - value) > 1e-10) { ret <- 1 }
  if (sign(ev$vectors[1,1]) != sign(vector[1])) { ev$vectors <- -ev$vectors }
  if (any(abs(ev$vectors[,1] - vector) > 1e-10)) { ret <- 1 }

  if (ret) { stop("Error") }

  0
}

g <- graph.famous("Zachary")
lc <- leading.eigenvector.community(g, callback=f)
lc
lc$modularity == modularity(g, lc$membership)
membership(lc)
length(lc)
sizes(lc)

## Check that the modularity matrix is correct

f <- function(membership, community, value, vector, multiplier, extra) {
  M <- sapply(1:length(vector), function(x) {
    v <- rep(0, length(vector))
    v[x] <- 1
    multiplier(v)
  })
  myc <- membership==community
  B <- A[myc,myc] - (deg[myc] %*% t(deg[myc]))/2/ec
  BG <- B-diag(rowSums(B))
  
  if (max(abs(M-BG)) > 1e-10) { stop("Error") }

  0
}

g <- graph.famous("Zachary")
A <- get.adjacency(g, sparse=FALSE)
ec <- ecount(g)
deg <- degree(g)
lc <- leading.eigenvector.community(g, callback=f)

## Stress-test

for (i in 1:100) {
  g <- erdos.renyi.game(20, sample(5:40, 1), type="gnm")
  lec1 <- leading.eigenvector.community(g)
  lec2 <- leading.eigenvector.community(g)
  if (length(membership(lec1)) != length(membership(lec2)) ||
      any(membership(lec1) != membership(lec2))) {
    print("Foobar!")
    break
  }
}

## Weighted

g <- graph.famous("Zachary")
lc1 <- leading.eigenvector.community(g)
lc2 <- leading.eigenvector.community(g, weights=rep(1, ecount(g)))

all(membership(lc1) == membership(lc2))
lc1$eigenvalues
lc2$eigenvalues

lc3 <- leading.eigenvector.community(g, weights=rep(1:4, length=ecount(g)))
membership(lc3)
