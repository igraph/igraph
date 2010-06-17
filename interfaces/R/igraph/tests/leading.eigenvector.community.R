
library(igraph)

g <- graph.famous("Zachary")
lc <- leading.eigenvector.community(g)
lc
lc$modularity == modularity(g, lc$membership)
membership(lc)
length(lc)
sizes(lc)
