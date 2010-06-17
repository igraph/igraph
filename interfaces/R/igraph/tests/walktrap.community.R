
library(igraph)

g <- graph.famous("Zachary")
set.seed(42)
wc <- walktrap.community(g)
wc
modularity(g, membership(wc)) == modularity(wc)
membership(wc)
modularity(wc)
length(wc)
sizes(wc)
d <- as.dendrogram(wc)
d
d[[1]]
d[[2]]
m2 <- cutat(wc, no=3)
modularity(g, m2) == wc$modularity[length(wc$modularity)-2]
