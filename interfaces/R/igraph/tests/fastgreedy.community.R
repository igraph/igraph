
library(igraph)

g <- graph.famous("Zachary")
fc <- fastgreedy.community(g)
fc
modularity(g, fc$membership) == max(fc$modularity)
membership(fc)
length(fc)
sizes(fc)
d <- as.dendrogram(fc)
d
d[[1]]
d[[2]]
m2 <- cutat(fc, no=3)
modularity(g, m2) == fc$modularity[length(fc$modularity)-2]
