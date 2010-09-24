
library(igraph)

set.seed(42)                            # Shouldn't be here....

ashs <- function(graph, as=TRUE) {
  mscale <- function(x) {
    if (sd(x)!=0) { x <- scale(x) }
    if (x[1] < 0) { x <- -x       }
    x
  }
  A <- get.adjacency(graph)
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




