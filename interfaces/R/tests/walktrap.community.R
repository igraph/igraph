
library(igraph)

g <- graph.famous("Zachary")
set.seed(42)
wc <- walktrap.community(g)
wc
abs(modularity(g, membership(wc)) - modularity(wc)) < 1e-7
membership(wc)
modularity(wc)
length(wc)
sizes(wc)
d <- as.dendrogram(wc)
d
d[[1]]
d[[2]]
m2 <- cutat(wc, no=3)
abs(modularity(g, m2) - wc$modularity[length(wc$modularity)-2]) < 1e-7
