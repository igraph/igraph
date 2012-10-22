
library(igraph)

g.1 <- graph( c(1,3,2,3,3,4,4,5) )
g.2 <- graph( c(2,1,3,1,4,1,5,1) )
g.3 <- graph( c(1,2,2,3,3,4,4,1,5,1) )
alpha.centrality(g.1)
alpha.centrality(g.2)
alpha.centrality(g.3,alpha=0.5)

alpha.centrality(g.1, sparse=FALSE)
alpha.centrality(g.2, sparse=FALSE)
alpha.centrality(g.3, alpha=0.5, sparse=FALSE)

##############################
## weighted version

set.seed(42)

star <- graph.star(10)
E(star)$weight <- runif(ecount(star))

alpha.centrality(star, sparse=TRUE)
alpha.centrality(star, sparse=FALSE)

alpha.centrality(star, weights="weight", sparse=TRUE)
alpha.centrality(star, weights="weight", sparse=FALSE)

alpha.centrality(star, weights=NA, sparse=TRUE)
alpha.centrality(star, weights=NA, sparse=FALSE)
