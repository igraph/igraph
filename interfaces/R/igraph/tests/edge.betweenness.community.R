
library(igraph)

g <- graph.famous("Zachary")
ebc <- edge.betweenness.community(g)
max(ebc$modularity) == modularity(g, ebc$membership)
ebc
membership(ebc)
length(ebc)
sizes(ebc)
d <- as.dendrogram(ebc)
d
d[[1]]
d[[2]]
