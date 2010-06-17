
library(igraph)

g <- graph.famous("Zachary")
oc <- optimal.community(g)
oc
membership(oc)
modularity(g, oc$membership) == oc$modularity

