
library(igraph)

g <- graph.famous("Zachary")
mc <- multilevel.community(g)
mc
membership(mc)
modularity(g, mc$membership) == max(mc$modularity)
length(mc)
sizes(mc)
