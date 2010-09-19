
library(igraph)

g.1 <- graph( c(1,3,2,3,3,4,4,5) )
g.2 <- graph( c(2,1,3,1,4,1,5,1) )
g.3 <- graph( c(1,2,2,3,3,4,4,1,5,1) )
alpha.centrality(g.1)
alpha.centrality(g.2)
alpha.centrality(g.3,alpha=0.5)

