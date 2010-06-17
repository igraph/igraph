
library(igraph)

g <- graph.famous("Zachary")
ebc <- edge.betweenness.community(g)
max(ebc$modularity) == modularity(g, ebc$membership)
ebc
membership(ebc)
