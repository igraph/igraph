
library(igraph)

g <- graph.famous("Zachary")
fc <- fastgreedy.community(g)
fc
abs(modularity(g, fc$membership) - max(fc$modularity)) < 1e-14
membership(fc)
length(fc)
sizes(fc)
d <- as.dendrogram(fc)
d
d[[1]]
d[[2]]
m2 <- cutat(fc, no=3)
abs(modularity(g, m2) - fc$modularity[length(fc$modularity)-2]) < 1e-14
