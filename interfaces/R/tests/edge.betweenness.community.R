
library(igraph)

g <- graph.famous("Zachary")
ebc <- edge.betweenness.community(g)
abs(max(ebc$modularity) - modularity(g, ebc$membership)) < 1e-14
ebc
membership(ebc)
length(ebc)
sizes(ebc)
d <- as.dendrogram(ebc)
d
d[[1]]
d[[2]]
m2 <- cutat(ebc, no=3)
abs(modularity(g, m2) - ebc$modularity[length(ebc$modularity)-2]) < 1e-14
