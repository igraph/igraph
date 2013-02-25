
library(igraph)

g <- graph.full(3)
path <- graph.ring(3, circular=F)
V(g)$color <- c(1,1,2)
V(path)$color <- c(1,2,1)
graph.count.subisomorphisms.vf2(g, path)
