
library(igraph)

g <- graph.famous("Zachary")
mc <- multilevel.community(g)
mc
membership(mc)
abs(modularity(g, mc$membership) - max(mc$modularity)) < 1e-14
length(mc)
sizes(mc)
